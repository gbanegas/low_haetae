#if defined(STM32F411E_DISCO)
#include <libopencm3/stm32/rcc.h>
#include <libopencm3/stm32/gpio.h>
#include <libopencm3/cm3/scb.h>
#include <libopencm3/cm3/dwt.h>
#endif

#include <stdlib.h>
#include <stdio.h>
#include <stdint.h>
#include <string.h>
#include <stdarg.h>

#include "hal.h"
#include "randombytes.h"
#include "m4.h"
#include "api.h"
#include "config.h"
#include "sign.h"
#include "keygen.h"
#include "stack_profile.h"

/* ---------------- Benchmark knobs ---------------- */
#define BENCH_CYCLES 1
#define BENCH_ITERS 10
/*
#ifndef BENCH_CYCLES
#define BENCH_CYCLES 1
#endif

#ifndef BENCH_ITERS
#define BENCH_ITERS 100
#endif
*/

/* STM32F411E-DISCO user LEDs:
 *   LD3: PD13 (orange)
 *   LD4: PD12 (green)
 *   LD5: PD14 (red)
 *   LD6: PD15 (blue)
 */
#if defined(STM32F411E_DISCO)
#define LED_PORT GPIOD
#define LED_KEYGEN GPIO12
#define LED_SIGN GPIO13
#define LED_VERIFY GPIO14
#define LED_ERR GPIO15

static void gpio_setup(void)
{
	rcc_periph_clock_enable(RCC_GPIOD);

	gpio_mode_setup(LED_PORT, GPIO_MODE_OUTPUT, GPIO_PUPD_NONE,
					LED_KEYGEN | LED_SIGN | LED_VERIFY | LED_ERR);
	gpio_set_output_options(LED_PORT, GPIO_OTYPE_PP, GPIO_OSPEED_2MHZ,
							LED_KEYGEN | LED_SIGN | LED_VERIFY | LED_ERR);
}

static void delay(volatile uint32_t n)
{
	while (n--)
	{
		__asm__ volatile("nop");
	}
}
#endif

static void hal_sendf(const char *fmt, ...)
{
	char buf[160];
	va_list ap;
	va_start(ap, fmt);
	vsnprintf(buf, sizeof(buf), fmt, ap);
	va_end(ap);
	hal_send_str(buf);
}

static void print_sizes(void)
{
	hal_sendf("N=%d K=%d L=%d", N, K, L);
	hal_sendf("sizeof(poly)=%zu", sizeof(poly));
	hal_sendf("sizeof(polyveck)=%zu", sizeof(polyveck));
	hal_sendf("sizeof(polyvecl)=%zu", sizeof(polyvecl));
	hal_sendf("sizeof(polyvecm)=%zu", sizeof(polyvecm));
}

/* ---------------- Cycles (DWT) ---------------- */
#if defined(STM32F411E_DISCO) && BENCH_CYCLES
#include <libopencm3/cm3/dwt.h>

static inline void cycles_init(void)
{
	dwt_enable_cycle_counter(); /* enables TRC + CYCCNT */
	// dwt_reset_cycle_counter();
}

static inline uint32_t cycles_now(void)
{
	return dwt_read_cycle_counter();
}
#endif

/* ---------------- Stats helpers ---------------- */
#if BENCH_CYCLES
static void insertion_sort_u32(uint32_t *a, size_t n)
{
	for (size_t i = 1; i < n; i++)
	{
		uint32_t x = a[i];
		size_t j = i;
		while (j > 0 && a[j - 1] > x)
		{
			a[j] = a[j - 1];
			j--;
		}
		a[j] = x;
	}
}

static uint32_t median_u32(const uint32_t *in, size_t n)
{
	/* n is small (100), O(n^2) sort is fine */
	static uint32_t tmp[BENCH_ITERS];
	if (n > BENCH_ITERS)
		n = BENCH_ITERS;

	for (size_t i = 0; i < n; i++)
		tmp[i] = in[i];
	insertion_sort_u32(tmp, n);

	if (n & 1)
	{
		return tmp[n / 2];
	}
	else
	{
		/* average of the two middle values */
		uint64_t a = tmp[(n / 2) - 1];
		uint64_t b = tmp[n / 2];
		return (uint32_t)((a + b) / 2);
	}
}

static void minmax_u32(const uint32_t *in, size_t n, uint32_t *mn, uint32_t *mx)
{
	uint32_t lo = 0xFFFFFFFFu, hi = 0;
	for (size_t i = 0; i < n; i++)
	{
		if (in[i] < lo)
			lo = in[i];
		if (in[i] > hi)
			hi = in[i];
	}
	*mn = lo;
	*mx = hi;
}
#endif

int main(void)
{
#if defined(STM32F411E_DISCO)
	gpio_setup();
	hal_setup(CLOCK_BENCHMARK);

	/* quick LED sweep */
	for (int i = 0; i < 2; i++)
	{
		gpio_toggle(LED_PORT, LED_KEYGEN);
		delay(200000);
		gpio_toggle(LED_PORT, LED_SIGN);
		delay(200000);
		gpio_toggle(LED_PORT, LED_VERIFY);
		delay(200000);
		gpio_toggle(LED_PORT, LED_ERR);
		delay(200000);
	}

#if BENCH_CYCLES
	cycles_init();
#endif
#endif

	hal_send_str("BOOT DONE!");
	hal_send_str("====== START ======");

	/* Put big buffers in .bss (not on stack) */
	static uint8_t pk[CRYPTO_PUBLICKEYBYTES];
	static uint8_t sk[CRYPTO_SECRETKEYBYTES];
	static uint8_t sig[CRYPTO_BYTES];
	static uint8_t m[32] = {0x1, 0x2, 0x3, 0x4};
	static size_t siglen = 0;
	const size_t mlen = sizeof(m);

	print_sizes();

	while (1)
	{
		hal_send_str("Starting!");

		/* ----------- Stack profiling pass (once) ----------- */
		sp_scope_t sc;

#if defined(STM32F411E_DISCO)
		gpio_set(LED_PORT, LED_KEYGEN);
#endif
		sc = sp_begin("crypto_sign_keypair_streamed_v2");
		crypto_sign_keypair_streamed_v2(pk, sk);
		sp_end(&sc);
#if defined(STM32F411E_DISCO)
		gpio_clear(LED_PORT, LED_KEYGEN);
#endif

#if defined(STM32F411E_DISCO)
		gpio_set(LED_PORT, LED_SIGN);
#endif
		sc = sp_begin("crypto_sign_signature");
		crypto_sign_signature(sig, &siglen, m, mlen, sk);
		sp_end(&sc);
#if defined(STM32F411E_DISCO)
		gpio_clear(LED_PORT, LED_SIGN);
#endif

#if defined(STM32F411E_DISCO)
		gpio_set(LED_PORT, LED_VERIFY);
#endif
		sc = sp_begin("crypto_sign_verify");
		int res = crypto_sign_verify(sig, siglen, m, mlen, pk);
		sp_end(&sc);
#if defined(STM32F411E_DISCO)
		gpio_clear(LED_PORT, LED_VERIFY);
#endif

		if (res != 0)
		{
			hal_send_str("VERIFY FAIL!");
#if defined(STM32F411E_DISCO)
			while (1)
			{
				gpio_toggle(LED_PORT, LED_ERR);
				delay(120000);
			}
#else
			while (1)
			{
			}
#endif
		}

		/* ----------- Cycles benchmark pass (median of 100) ----------- */
#if defined(STM32F411E_DISCO) && BENCH_CYCLES
		static uint32_t cyc_key[BENCH_ITERS];
		static uint32_t cyc_sign[BENCH_ITERS];
		static uint32_t cyc_vfy[BENCH_ITERS];

		for (size_t i = 0; i < BENCH_ITERS; i++)
		{
			uint32_t t0, t1;

			t0 = cycles_now();
			crypto_sign_keypair_streamed_v2(pk, sk);
			t1 = cycles_now();
			cyc_key[i] = (uint32_t)(t1 - t0);

			for (int i = 0; i < 2; i++)
			{
				gpio_toggle(LED_PORT, LED_KEYGEN);
				delay(200000);
				gpio_toggle(LED_PORT, LED_SIGN);
				delay(200000);
				gpio_toggle(LED_PORT, LED_VERIFY);
				delay(200000);
				
			}

			t0 = cycles_now();
			crypto_sign_signature(sig, &siglen, m, mlen, sk);
			t1 = cycles_now();
			cyc_sign[i] = (uint32_t)(t1 - t0);

			for (int i = 0; i < 2; i++)
			{
				gpio_toggle(LED_PORT, LED_KEYGEN);
				delay(200000);
				gpio_toggle(LED_PORT, LED_SIGN);
				delay(200000);
				gpio_toggle(LED_PORT, LED_VERIFY);
				delay(200000);
				
			}

			t0 = cycles_now();
			res = crypto_sign_verify(sig, siglen, m, mlen, pk);
			t1 = cycles_now();
			cyc_vfy[i] = (uint32_t)(t1 - t0);

			if (res != 0)
			{
				hal_sendf("VERIFY FAIL during bench at i=%u", (unsigned)i);
				while (1)
				{
					gpio_toggle(LED_PORT, LED_ERR);
					delay(120000);
				}
			}
		}

		uint32_t med_k = median_u32(cyc_key, BENCH_ITERS);
		uint32_t med_s = median_u32(cyc_sign, BENCH_ITERS);
		uint32_t med_v = median_u32(cyc_vfy, BENCH_ITERS);

		uint32_t min_k, max_k, min_s, max_s, min_v, max_v;
		minmax_u32(cyc_key, BENCH_ITERS, &min_k, &max_k);
		minmax_u32(cyc_sign, BENCH_ITERS, &min_s, &max_s);
		minmax_u32(cyc_vfy, BENCH_ITERS, &min_v, &max_v);

		hal_sendf("[CY] iters=%u", (unsigned)BENCH_ITERS);
		hal_sendf("[CY] keygen  median=%u  min=%u  max=%u", (unsigned)med_k, (unsigned)min_k, (unsigned)max_k);
		hal_sendf("[CY] sign    median=%u  min=%u  max=%u", (unsigned)med_s, (unsigned)min_s, (unsigned)max_s);
		hal_sendf("[CY] verify  median=%u  min=%u  max=%u", (unsigned)med_v, (unsigned)min_v, (unsigned)max_v);
#endif

		/* idle marker */
#if defined(STM32F411E_DISCO)
		gpio_toggle(LED_PORT, LED_ERR);
		delay(800000);
#endif
	}
}