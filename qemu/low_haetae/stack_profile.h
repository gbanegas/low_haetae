#ifndef STACK_PROFILE_H
#define STACK_PROFILE_H

#include <stdint.h>
#include <stddef.h>
#include "hal.h"
#define STACK_PROFILE 1


#define STACK_PROFILE_USE_HAL 1


extern uint8_t __stack_bottom__;
extern uint8_t __stack_top__;

typedef struct {
    const char *tag;
    uintptr_t sp_entry;
    uintptr_t bottom;
    uintptr_t top;
    uint32_t base_used;
    uint32_t peak_used;
} sp_scope_t;

#if STACK_PROFILE

static inline __attribute__((always_inline)) uintptr_t sp_read(void) {
    uintptr_t sp;
    __asm__ volatile ("mov %0, sp" : "=r"(sp));
    return sp;
}

static inline __attribute__((always_inline)) uint32_t sp_used_from(uintptr_t sp) {
    uintptr_t bottom = (uintptr_t)&__stack_bottom__;
    uintptr_t top    = (uintptr_t)&__stack_top__;
    if (sp < bottom) sp = bottom;
    if (sp > top)    sp = top;
    return (uint32_t)(top - sp);
}

static inline __attribute__((always_inline)) void sp_paint(uintptr_t bottom, uintptr_t sp_entry) {
    /* paint 32-bit words for speed */
    uint32_t *p = (uint32_t*)bottom;
    uint32_t *e = (uint32_t*)sp_entry;
    while (p < e) *p++ = 0xA5A5A5A5u;
}

/* returns deepest-used SP (lowest address) */
static inline __attribute__((always_inline)) uintptr_t sp_scan_deepest(uintptr_t bottom, uintptr_t sp_entry) {
    uint32_t *p = (uint32_t*)bottom;
    uint32_t *e = (uint32_t*)sp_entry;

    while (p < e && *p == 0xA5A5A5A5u) p++;

    /* byte-precision refinement inside the first dirty word */
    uintptr_t addr = (uintptr_t)p;
    if (addr >= sp_entry) return sp_entry;

    uint8_t *b = (uint8_t*)addr;
    for (int i = 0; i < 4; i++) {
        if (b[i] != 0xA5) return (uintptr_t)(b + i);
    }
    return addr;
}

void sp_print(const sp_scope_t *s); /* implemented in .c */

static inline __attribute__((always_inline)) sp_scope_t sp_begin(const char *tag) {
    sp_scope_t s;
    s.tag     = tag;
    s.bottom  = (uintptr_t)&__stack_bottom__;
    s.top     = (uintptr_t)&__stack_top__;
    s.sp_entry= sp_read();
    s.base_used = sp_used_from(s.sp_entry);
    s.peak_used = s.base_used;

    sp_paint(s.bottom, s.sp_entry);
    return s;
}

static inline __attribute__((always_inline)) void sp_end(sp_scope_t *s) {
    /* scan BEFORE printing to avoid profiler affecting result */
    uintptr_t deepest = sp_scan_deepest(s->bottom, s->sp_entry);
    uint32_t peak = sp_used_from(deepest);
    if (peak > s->peak_used) s->peak_used = peak;

    sp_print(s);
}

#else

static inline sp_scope_t sp_begin(const char *tag) { (void)tag; sp_scope_t s={0}; return s; }
static inline void sp_end(sp_scope_t *s) { (void)s; }

#endif

#endif
