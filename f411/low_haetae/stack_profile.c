#include "stack_profile.h"
#include <stdint.h>
#include <stddef.h>

#if STACK_PROFILE

static inline char *app_str(char *p, const char *s) {
    while (*s) *p++ = *s++;
    return p;
}
static inline char *app_u32(char *p, uint32_t v) {
    char tmp[10];
    int n = 0;
    if (v == 0) { *p++ = '0'; return p; }
    while (v && n < 10) { tmp[n++] = (char)('0' + (v % 10)); v /= 10; }
    while (n--) *p++ = tmp[n];
    return p;
}

static inline uintptr_t read_sp(void)
{
    uintptr_t sp;
    __asm__ volatile ("mov %0, sp" : "=r"(sp));
    return sp;
}

#ifndef SP_GUARD_BYTES
#define SP_GUARD_BYTES 256u
#endif

static void paint_a5(uintptr_t lo, uintptr_t hi)
{
    volatile uint8_t *p = (volatile uint8_t *)lo;
    volatile uint8_t *e = (volatile uint8_t *)hi;
    while (p < e) *p++ = 0xA5;
}

static uintptr_t scan_first_non_a5(uintptr_t lo, uintptr_t hi)
{
    volatile uint8_t *p = (volatile uint8_t *)lo;
    volatile uint8_t *e = (volatile uint8_t *)hi;
    while (p < e) {
        if (*p != 0xA5) return (uintptr_t)p;
        p++;
    }
    return hi; /* untouched */
}

static void sp_print(const sp_scope_t *s)
{
#if STACK_PROFILE_USE_HAL
    static char buf[140];
    char *p = buf;
    uint32_t delta = (s->peak >= s->base) ? (s->peak - s->base) : 0u;

    p = app_str(p, "[SP] ");
    p = app_str(p, s->name ? s->name : "?");
    p = app_str(p, " base=");
    p = app_u32(p, s->base);
    p = app_str(p, " peak=");
    p = app_u32(p, s->peak);
    p = app_str(p, " (+");
    p = app_u32(p, delta);
    p = app_str(p, ") / ");
    p = app_u32(p, s->span);
    *p = 0;
    hal_send_str(buf);
#else
    (void)s;
#endif
}

sp_scope_t sp_begin(const char *name)
{
    sp_scope_t s;
    s.name = name;

    s.span = sp_total_span();
    s.base = sp_used_now();
    s.peak = s.base;

    uintptr_t lo = (uintptr_t)&end;
    uintptr_t sp = read_sp();
    uintptr_t hi = (sp > SP_GUARD_BYTES) ? (sp - SP_GUARD_BYTES) : sp;

    s.paint_lo = lo;
    s.paint_hi = hi;

    if (hi > lo) paint_a5(lo, hi);

    return s;
}

void sp_sample(sp_scope_t *s)
{
    /* boundary sampling helper: only useful if called INSIDE long scopes */
    uint32_t now = sp_used_now();
    if (now > s->peak) s->peak = now;
}

void sp_end(sp_scope_t *s)
{
    uintptr_t lo  = s->paint_lo;
    uintptr_t hi  = s->paint_hi;

    if (hi > lo) {
        uintptr_t first = scan_first_non_a5(lo, hi);

        /* If untouched, don't “charge” the guard as peak */
        if (first != hi) {
            uintptr_t top = (uintptr_t)&_stack;
            uint32_t used = (top > first) ? (uint32_t)(top - first) : 0u;
            if (used > s->peak) s->peak = used;
        }
    }

    sp_print(s);
}

#endif