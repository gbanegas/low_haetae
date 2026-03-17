\
#include <stdint.h>

/* Linker-provided symbols */
extern uint32_t _estack;

extern uint32_t _sidata; /* start of init values in FLASH */
extern uint32_t _sdata;  /* start of .data in RAM */
extern uint32_t _edata;  /* end of .data in RAM */
extern uint32_t _sbss;   /* start of .bss in RAM */
extern uint32_t _ebss;   /* end of .bss in RAM */

int main(void);

/* Default handlers */
void reset_handler(void);
void default_handler(void);

void __attribute__((weak, alias("default_handler"))) nmi_handler(void);
void __attribute__((weak, alias("default_handler"))) hard_fault_handler(void);
void __attribute__((weak, alias("default_handler"))) mem_manage_handler(void);
void __attribute__((weak, alias("default_handler"))) bus_fault_handler(void);
void __attribute__((weak, alias("default_handler"))) usage_fault_handler(void);
void __attribute__((weak, alias("default_handler"))) svc_handler(void);
void __attribute__((weak, alias("default_handler"))) debug_monitor_handler(void);
void __attribute__((weak, alias("default_handler"))) pendsv_handler(void);
void __attribute__((weak, alias("default_handler"))) sys_tick_handler(void);

/* Minimal vector table in the section libopencm3 commonly uses. */
__attribute__((section(".vectors")))
void (*const vector_table[])(void) = {
    (void (*)(void))(&_estack), /* Initial stack pointer */
    reset_handler,
    nmi_handler,
    hard_fault_handler,
    mem_manage_handler,
    bus_fault_handler,
    usage_fault_handler,
    0, 0, 0, 0,                /* Reserved */
    svc_handler,
    debug_monitor_handler,
    0,                         /* Reserved */
    pendsv_handler,
    sys_tick_handler,
    /* Device IRQs not needed for this blinky. */
};

void reset_handler(void)
{
    /* Copy .data from FLASH to RAM */
    uint32_t *src = &_sidata;
    uint32_t *dst = &_sdata;
    while (dst < &_edata) {
        *dst++ = *src++;
    }

    /* Zero .bss */
    dst = &_sbss;
    while (dst < &_ebss) {
        *dst++ = 0;
    }

    (void)main();

    while (1) {
        __asm__ volatile ("wfi");
    }
}

void default_handler(void)
{
    while (1) {
        __asm__ volatile ("wfi");
    }
}
