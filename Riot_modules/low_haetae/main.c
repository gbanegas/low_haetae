/*
 * RIOT benchmark harness for HAETAE: keygen, sign, verify
 *
 * Adds per-function peak stack usage.
 *
 * IMPORTANT (Cortex-M + RIOT mpu_stack_guard):
 *   RIOT protects a small region at the bottom of each thread stack with the MPU
 *   to catch overflows. Writing to the very first bytes of the stack can trigger
 *   a MemManage fault.
 *
 *   This file therefore SKIPS the first STACK_GUARD_SKIP_BYTES bytes when
 *   painting/scanning.
 */

#include <stdio.h>
#include <stdint.h>
#include <stdbool.h>
#include <stddef.h>
#include <string.h>
#include <inttypes.h>
#include <stdarg.h>

#include "thread.h"
#include "xtimer.h"
#include "msg.h"
#include "mutex.h"
#include "clk.h"
#include "api.h"
#include "keygen.h"
#include "sign.h"
#include "verify.h"

/* ---------------- configuration ---------------- */
#ifndef BENCH_KEYGEN_ITERS
#define BENCH_KEYGEN_ITERS  3
#endif
#ifndef BENCH_SIGN_ITERS
#define BENCH_SIGN_ITERS    10
#endif
#ifndef BENCH_VERIFY_ITERS
#define BENCH_VERIFY_ITERS  10
#endif

#ifndef STACK_PAINT_MARGIN_BYTES
/* Leave some bytes unpainted below the current SP to avoid painting into this
 * measurement routine's own transient stack usage. */
#define STACK_PAINT_MARGIN_BYTES  128u
#endif

#ifndef STACK_GUARD_SKIP_BYTES
/* Conservative default: RIOT's MPU stack guard commonly protects at least 32B.
 * 256B keeps us well away from MPU-aligned guard regions on Cortex-M4.
 * You can override with CFLAGS+='-DSTACK_GUARD_SKIP_BYTES=...'
 */
#define STACK_GUARD_SKIP_BYTES  256u
#endif

#define STACK_MARK  0xA5u

/* ---------------- print locking ---------------- */
static mutex_t print_lock = MUTEX_INIT;
static inline void safe_printf(const char *fmt, ...)
{
    va_list ap;
    mutex_lock(&print_lock);
    va_start(ap, fmt);
    vprintf(fmt, ap);
    va_end(ap);
    fflush(stdout);
    mutex_unlock(&print_lock);
}

/* ---------------- stack peak measurement ---------------- */
static inline uintptr_t _get_sp(void)
{
    uintptr_t sp;
#if defined(__riscv)
    __asm__ volatile ("mv %0, sp" : "=r"(sp));
#elif defined(__arm__)
    __asm__ volatile ("mov %0, sp" : "=r"(sp));
#elif defined(__xtensa__)
    __asm__ volatile ("mov %0, a1" : "=r"(sp));
#else
    /* Portable fallback: take address of a local and adjust for frame */
    volatile uintptr_t dummy = 0;
    sp = (uintptr_t)&dummy;
#endif
    return sp;
}

typedef struct {
    uint8_t   *stack_start;   /* lowest address */
    uint8_t   *stack_end;     /* highest address */
    uint8_t   *paint_start;   /* start of painted region */
    uint8_t   *paint_end;     /* end of painted region */
    uintptr_t  sp0;
} stack_ctx_t;

static inline void _stack_ctx_begin(stack_ctx_t *ctx, const thread_t *thr)
{
    uint8_t *ss = (uint8_t *)thread_get_stackstart(thr);
    size_t   sz = thread_get_stacksize(thr);
    uint8_t *se = ss + sz;

    ctx->stack_start = ss;
    ctx->stack_end   = se;

    /* snapshot SP before running the measured function */
    ctx->sp0 = _get_sp();

    /* Skip bottom guard area */
    uintptr_t ps = (uintptr_t)ss + (uintptr_t)STACK_GUARD_SKIP_BYTES;
    if (ps > (uintptr_t)se) {
        ps = (uintptr_t)se;
    }
    ctx->paint_start = (uint8_t *)ps;

    /* paint end: a bit below the current SP */
    uintptr_t pe = ctx->sp0;
    if (pe > (uintptr_t)ctx->paint_start + STACK_PAINT_MARGIN_BYTES) {
        pe -= STACK_PAINT_MARGIN_BYTES;
    }
    else {
        pe = (uintptr_t)ctx->paint_start;
    }

    if (pe > (uintptr_t)se) {
        pe = (uintptr_t)se;
    }

    ctx->paint_end = (uint8_t *)pe;

    /* Paint [paint_start, paint_end) */
    for (uint8_t *p = ctx->paint_start; p < ctx->paint_end; p++) {
        *p = (uint8_t)STACK_MARK;
    }
}

static inline void _stack_ctx_end(const stack_ctx_t *ctx,
                                 size_t *peak_used_total_B,
                                 size_t *delta_below_sp0_B)
{
    uint8_t *p = ctx->paint_start;
    while (p < ctx->paint_end && *p == (uint8_t)STACK_MARK) {
        p++;
    }

    /* If nothing overwritten in painted region, we only know that usage didn't
     * dip below (sp0 - margin). Use baseline estimate from sp0. */
    if (p == ctx->paint_end) {
        const uint8_t *sp0p = (const uint8_t *)ctx->sp0;
        if (sp0p >= ctx->stack_start && sp0p <= ctx->stack_end) {
            *peak_used_total_B = (size_t)(ctx->stack_end - sp0p);
        }
        else {
            *peak_used_total_B = 0;
        }
        *delta_below_sp0_B = 0;
        return;
    }

    /* first non-marker byte is the deepest point reached (lowest addr) */
    *peak_used_total_B = (size_t)(ctx->stack_end - p);

    const uint8_t *sp0p = (const uint8_t *)ctx->sp0;
    if (sp0p >= p) {
        *delta_below_sp0_B = (size_t)(sp0p - p);
    }
    else {
        *delta_below_sp0_B = 0;
    }
}

/* ---------------- bench bookkeeping ---------------- */
typedef struct {
    uint32_t min_us, max_us;
    uint64_t sum_us;
#ifdef CLOCK_CORECLOCK
    uint32_t min_ticks, max_ticks;
    uint64_t sum_ticks;
#endif
    size_t   stack_peak_max_B;
    size_t   stack_delta_max_B;
    uint32_t iters;
} bench_stat_t;

static inline void _bench_init(bench_stat_t *b)
{
    memset(b, 0, sizeof(*b));
    b->min_us = UINT32_MAX;
#ifdef CLOCK_CORECLOCK
    b->min_ticks = UINT32_MAX;
#endif
}

static inline void _bench_update(bench_stat_t *b, uint32_t dt_us,
                                 size_t peak_stack_B, size_t delta_stack_B)
{
    if (dt_us < b->min_us) b->min_us = dt_us;
    if (dt_us > b->max_us) b->max_us = dt_us;
    b->sum_us += dt_us;
#ifdef CLOCK_CORECLOCK
    uint32_t dt_ticks = (uint32_t)(((uint64_t)dt_us * (uint64_t)CLOCK_CORECLOCK) / (uint64_t)US_PER_SEC);
    if (dt_ticks < b->min_ticks) b->min_ticks = dt_ticks;
    if (dt_ticks > b->max_ticks) b->max_ticks = dt_ticks;
    b->sum_ticks += dt_ticks;
#endif
    if (peak_stack_B > b->stack_peak_max_B) b->stack_peak_max_B = peak_stack_B;
    if (delta_stack_B > b->stack_delta_max_B) b->stack_delta_max_B = delta_stack_B;
    b->iters++;
}

static inline uint32_t _bench_avg_us(const bench_stat_t *b)
{
    return (b->iters == 0) ? 0u : (uint32_t)(b->sum_us / b->iters);
}
#ifdef CLOCK_CORECLOCK
static inline uint32_t _bench_avg_ticks(const bench_stat_t *b)
{
    return (b->iters == 0) ? 0u : (uint32_t)(b->sum_ticks / b->iters);
}
#endif

/* ---------------- worker thread ---------------- */
#ifndef HAETAE_THREAD_STACKSIZE
#define HAETAE_THREAD_STACKSIZE (64u * 1024u)
#endif

static uint8_t pk[CRYPTO_PUBLICKEYBYTES];
static uint8_t sk[CRYPTO_SECRETKEYBYTES];

static bench_stat_t bench_keygen;
static bench_stat_t bench_sign;
static bench_stat_t bench_verify;

static size_t last_siglen = 0;
static int verify_fail = 0;

static char haetae_stack[HAETAE_THREAD_STACKSIZE] __attribute__((aligned(8)));
static kernel_pid_t main_pid;

static void *haetae_worker(void *arg)
{
    (void)arg;

    safe_printf("[haetae] worker start (stack=%u, guard_skip=%u, margin=%u)\n",
                (unsigned)sizeof(haetae_stack),
                (unsigned)STACK_GUARD_SKIP_BYTES,
                (unsigned)STACK_PAINT_MARGIN_BYTES);

    _bench_init(&bench_keygen);
    _bench_init(&bench_sign);
    _bench_init(&bench_verify);

    static const uint8_t msgbuf[32] = {
        0x48,0x41,0x45,0x54,0x41,0x45,0x2d,0x52,
        0x49,0x4f,0x54,0x2d,0x42,0x45,0x4e,0x43,
        0x48,0x2d,0x4d,0x53,0x47,0x2d,0x30,0x31,
        0x32,0x33,0x34,0x35,0x36,0x37,0x38,0x39
    };

    uint8_t sig[CRYPTO_BYTES];
    size_t siglen = 0;

    const thread_t *self = thread_get_active();

    /* ---- keygen ---- */
    for (unsigned i = 0; i < (unsigned)BENCH_KEYGEN_ITERS; i++) {
        stack_ctx_t sctx;
        _stack_ctx_begin(&sctx, self);
        uint32_t t0 = xtimer_now_usec();
        (void)crypto_sign_keypair_low_stack(pk, sk);
        uint32_t t1 = xtimer_now_usec();
        size_t peakB = 0, deltaB = 0;
        _stack_ctx_end(&sctx, &peakB, &deltaB);
        _bench_update(&bench_keygen, (uint32_t)(t1 - t0), peakB, deltaB);
        thread_yield();
    }

    /* ---- sign ---- */
    for (unsigned i = 0; i < (unsigned)BENCH_SIGN_ITERS; i++) {
        stack_ctx_t sctx;
        _stack_ctx_begin(&sctx, self);
        uint32_t t0 = xtimer_now_usec();
        int rc = crypto_sign_signature(sig, &siglen, msgbuf, sizeof(msgbuf), sk);
        uint32_t t1 = xtimer_now_usec();
        size_t peakB = 0, deltaB = 0;
        _stack_ctx_end(&sctx, &peakB, &deltaB);
        _bench_update(&bench_sign, (uint32_t)(t1 - t0), peakB, deltaB);

        if (rc != 0) {
            verify_fail++;
        }
        else if (crypto_sign_verify(sig, siglen, msgbuf, sizeof(msgbuf), pk) != 0) {
            verify_fail++;
        }
        thread_yield();
    }
    last_siglen = siglen;

    if (last_siglen == 0) {
        (void)crypto_sign_signature(sig, &siglen, msgbuf, sizeof(msgbuf), sk);
        last_siglen = siglen;
    }

    /* ---- verify ---- */
    for (unsigned i = 0; i < (unsigned)BENCH_VERIFY_ITERS; i++) {
        stack_ctx_t sctx;
        _stack_ctx_begin(&sctx, self);
        uint32_t t0 = xtimer_now_usec();
        int ok = crypto_sign_verify(sig, last_siglen, msgbuf, sizeof(msgbuf), pk);
        uint32_t t1 = xtimer_now_usec();
        size_t peakB = 0, deltaB = 0;
        _stack_ctx_end(&sctx, &peakB, &deltaB);
        _bench_update(&bench_verify, (uint32_t)(t1 - t0), peakB, deltaB);
        if (ok != 0) {
            verify_fail++;
        }
        thread_yield();
    }

    safe_printf("[haetae] worker done\n");

    msg_t m;
    m.type = 0xBEEF;
    msg_send(&m, main_pid);

    return NULL;
}

/* ---------------- main ---------------- */
static msg_t main_queue[8];

int main(void)
{
    puts("main(): This is RIOT! (HAETAE bench)");

    msg_init_queue(main_queue, (unsigned) (sizeof(main_queue) / sizeof(main_queue[0])));
    main_pid = thread_getpid();

    kernel_pid_t pid = thread_create(
        haetae_stack, sizeof(haetae_stack),
        THREAD_PRIORITY_MAIN - 1,
        THREAD_CREATE_STACKTEST,
        haetae_worker, NULL, "haetae");

    if (pid <= KERNEL_PID_UNDEF) {
        puts("thread_create failed!");
        return 1;
    }

    msg_t m;
    msg_receive(&m);

    /* Print one JSON line */
    printf("{\"bench\":{");

    printf("\"keygen\":{\"iters\":%" PRIu32 ",\"min_us\":%" PRIu32 ",\"avg_us\":%" PRIu32 ",\"max_us\":%" PRIu32 ","
           "\"stack_peak_B\":%lu,\"stack_delta_B\":%lu"
#ifdef CLOCK_CORECLOCK
           ",\"min_ticks\":%" PRIu32 ",\"avg_ticks\":%" PRIu32 ",\"max_ticks\":%" PRIu32
#endif
           "},",
           bench_keygen.iters,
           bench_keygen.min_us, _bench_avg_us(&bench_keygen), bench_keygen.max_us,
           (unsigned long)bench_keygen.stack_peak_max_B,
           (unsigned long)bench_keygen.stack_delta_max_B
#ifdef CLOCK_CORECLOCK
           , bench_keygen.min_ticks, _bench_avg_ticks(&bench_keygen), bench_keygen.max_ticks
#endif
    );

    printf("\"sign\":{\"iters\":%" PRIu32 ",\"min_us\":%" PRIu32 ",\"avg_us\":%" PRIu32 ",\"max_us\":%" PRIu32 ","
           "\"siglen\":%lu,\"stack_peak_B\":%lu,\"stack_delta_B\":%lu"
#ifdef CLOCK_CORECLOCK
           ",\"min_ticks\":%" PRIu32 ",\"avg_ticks\":%" PRIu32 ",\"max_ticks\":%" PRIu32
#endif
           "},",
           bench_sign.iters,
           bench_sign.min_us, _bench_avg_us(&bench_sign), bench_sign.max_us,
           (unsigned long)last_siglen,
           (unsigned long)bench_sign.stack_peak_max_B,
           (unsigned long)bench_sign.stack_delta_max_B
#ifdef CLOCK_CORECLOCK
           , bench_sign.min_ticks, _bench_avg_ticks(&bench_sign), bench_sign.max_ticks
#endif
    );

    printf("\"verify\":{\"iters\":%" PRIu32 ",\"min_us\":%" PRIu32 ",\"avg_us\":%" PRIu32 ",\"max_us\":%" PRIu32 ","
           "\"fail\":%d,\"stack_peak_B\":%lu,\"stack_delta_B\":%lu"
#ifdef CLOCK_CORECLOCK
           ",\"min_ticks\":%" PRIu32 ",\"avg_ticks\":%" PRIu32 ",\"max_ticks\":%" PRIu32
#endif
           "}",
           bench_verify.iters,
           bench_verify.min_us, _bench_avg_us(&bench_verify), bench_verify.max_us,
           verify_fail,
           (unsigned long)bench_verify.stack_peak_max_B,
           (unsigned long)bench_verify.stack_delta_max_B
#ifdef CLOCK_CORECLOCK
           , bench_verify.min_ticks, _bench_avg_ticks(&bench_verify), bench_verify.max_ticks
#endif
    );

    printf("}");
#ifdef CLOCK_CORECLOCK
    printf(",\"coreclock_hz\":%" PRIu32, (uint32_t)CLOCK_CORECLOCK);
#endif
    printf(",\"thread_stack_B\":%u,\"guard_skip_B\":%u,\"margin_B\":%u}\n",
           (unsigned)sizeof(haetae_stack),
           (unsigned)STACK_GUARD_SKIP_BYTES,
           (unsigned)STACK_PAINT_MARGIN_BYTES);
    fflush(stdout);

    while (1) {
        xtimer_sleep(1);
    }
}
