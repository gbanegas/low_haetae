#include <libopencm3/stm32/rcc.h>
#include <libopencm3/stm32/gpio.h>

void hard_fault_handler(void)
{
    /* Blue LED LD6 on STM32F411E-DISCO is PD15. */
    rcc_periph_clock_enable(RCC_GPIOD);
    gpio_mode_setup(GPIOD, GPIO_MODE_OUTPUT, GPIO_PUPD_NONE, GPIO15);
    gpio_set_output_options(GPIOD, GPIO_OTYPE_PP, GPIO_OSPEED_2MHZ, GPIO15);

    while (1) {
        gpio_toggle(GPIOD, GPIO15);
        for (volatile unsigned i = 0; i < 120000; i++) __asm__ volatile ("nop");
    }
}
