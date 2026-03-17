#include "stack_profile.h"

#if STACK_PROFILE && STACK_PROFILE_USE_HAL
#include "hal.h"
#endif

#if STACK_PROFILE

static inline char *app_str(char *p, const char *s) {
    while (*s) *p++ = *s++;
    return p;
}
static inline char *app_u32(char *p, uint32_t v) {
    char tmp[10];
    int n = 0;
    if (v == 0) { *p++='0'; return p; }
    while (v && n < 10) { tmp[n++] = (char)('0' + (v % 10)); v /= 10; }
    while (n--) *p++ = tmp[n];
    return p;
}

void sp_print(const sp_scope_t *s) {
#if STACK_PROFILE_USE_HAL
    static char buf[128];
    char *p = buf;
    uint32_t delta = (s->peak_used >= s->base_used) ? (s->peak_used - s->base_used) : 0;
    uint32_t span  = (uint32_t)(s->top - s->bottom);

    p = app_str(p, "[SP] ");
    p = app_str(p, s->tag ? s->tag : "?");
    p = app_str(p, " base=");
    p = app_u32(p, s->base_used);
    p = app_str(p, " peak=");
    p = app_u32(p, s->peak_used);
    p = app_str(p, " (+");
    p = app_u32(p, delta);
    p = app_str(p, ") / ");
    p = app_u32(p, span);
    *p = 0;

    hal_send_str(buf);
#else
    (void)s;
#endif
}

#endif
