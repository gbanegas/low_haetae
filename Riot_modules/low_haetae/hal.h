#ifndef HAETAE_HAL_H
#define HAETAE_HAL_H

/* Minimal HAL used by the demo + stack profiler.
 * On RIOT, this is a thin wrapper around stdio.
 */
void hal_send_str(const char *s);

#endif /* HAETAE_HAL_H */
