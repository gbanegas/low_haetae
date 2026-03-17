#ifndef STACK_PROFILE_H
#define STACK_PROFILE_H

#include <stdint.h>
#include <stddef.h>
#include "hal.h"

#define STACK_PROFILE 1
#define STACK_PROFILE_USE_HAL 1

typedef struct {
    const char *name;

    /* Painted window [paint_lo, paint_hi) at begin */
    uintptr_t paint_lo;
    uintptr_t paint_hi;

    /* Total stack span = &_stack - &end */
    uint32_t span;

    /* Measurements in bytes "used from top" */
    uint32_t base;   /* used at begin */
    uint32_t peak;   /* max used during scope */
} sp_scope_t;

#if STACK_PROFILE

/* Provided by linker script */
extern char end;     /* end of .bss / start of heap */
extern char _stack;  /* initial stack pointer (top of stack) */

/* total available stack bytes */
static inline uint32_t sp_total_span(void)
{
    uintptr_t lo = (uintptr_t)&end;
    uintptr_t top = (uintptr_t)&_stack;
    return (top > lo) ? (uint32_t)(top - lo) : 0u;
}

/* used stack right now (bytes from &_stack down to current SP) */
static inline uint32_t sp_used_now(void)
{
    uint32_t span = sp_total_span();
    uint32_t free = (uint32_t)hal_get_stack_size(); /* SP - heap_start */
    if (free > span) free = span;
    return span - free;
}

sp_scope_t sp_begin(const char *name);
void sp_sample(sp_scope_t *s);   /* optional mid-scope sampling */
void sp_end(sp_scope_t *s);

#else

static inline sp_scope_t sp_begin(const char *name) { (void)name; sp_scope_t s = {0}; return s; }
static inline void sp_sample(sp_scope_t *s) { (void)s; }
static inline void sp_end(sp_scope_t *s) { (void)s; }

#endif

#endif