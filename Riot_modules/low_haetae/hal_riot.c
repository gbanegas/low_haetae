#include "hal.h"

#include <stdio.h>
#include <string.h>

void hal_send_str(const char *s)
{
    if (!s) {
        return;
    }

    /* RIOT's puts() appends a newline. Avoid double newlines when the caller
       already provides one. */
    size_t n = strlen(s);
    if (n && s[n - 1] == '\n') {
        fputs(s, stdout);
    }
    else {
        puts(s);
    }
}
