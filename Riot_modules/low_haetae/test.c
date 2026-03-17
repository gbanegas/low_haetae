
/*
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
#include "stack_profile.h"

static void fill_message(uint8_t *m, size_t mlen)
{
	for (size_t i = 0; i < mlen; i++)
	{
		m[i] = (uint8_t)i;
	}
}

static void print_bytes_dec(const char *label, const uint8_t *buf, size_t len)
{
	char out[256];
	size_t pos = 0;

	// Optional label
	if (label)
	{
		pos += snprintf(out + pos, sizeof(out) - pos, "%s", label);
	}

	for (size_t i = 0; i < len; i++)
	{
		// If we’re close to the end of the buffer, flush
		if (pos > sizeof(out) - 8)
		{
			out[pos] = '\0';
			hal_send_str(out);
			pos = 0;
		}
		pos += snprintf(out + pos, sizeof(out) - pos, "%u,", buf[i]);
	}

	// Add newline and flush remaining
	if (pos > sizeof(out) - 4)
	{
		out[pos] = '\0';
		hal_send_str(out);
		pos = 0;
	}
	pos += snprintf(out + pos, sizeof(out) - pos, "\n");
	out[pos] = '\0';
	hal_send_str(out);
}

static void hal_sendf(const char *fmt, ...)
{
	char buf[128]; // increase if you want longer lines
	va_list ap;
	va_start(ap, fmt);
	vsnprintf(buf, sizeof(buf), fmt, ap);
	va_end(ap);

	hal_send_str(buf); // hal_send_str adds '\n' already
}

void print_sizes(void)
{
	hal_sendf("N=%d K=%d L=%d", N, K, L);
	hal_sendf("sizeof(poly)=%zu", sizeof(poly));
	hal_sendf("sizeof(polyveck)=%zu", sizeof(polyveck));
	hal_sendf("sizeof(polyvecl)=%zu", sizeof(polyvecl));
	hal_sendf("sizeof(polyvecm)=%zu", sizeof(polyvecm));
}

extern uint8_t __stack_bottom__, __stack_top__;

int main(void)
{

	sp_scope_t sc = sp_begin("print");

	hal_send_str("====== START ======");

	sp_end(&sc);
	print_sizes();

	hal_sendf("stack span = %lu", (unsigned long)(&__stack_top__ - &__stack_bottom__));

	uint8_t pk[CRYPTO_PUBLICKEYBYTES];
	uint8_t sk[CRYPTO_SECRETKEYBYTES];

	// Small test message
	uint8_t m[32] = {0x1, 0x2, 0x3, 0x4};
	size_t mlen = sizeof m;

	// Signature buffer
	uint8_t sig[CRYPTO_BYTES];
	size_t siglen = 0;

	while (1)
	{
		hal_send_str("Starting the KeyPair !\n");

		sc = sp_begin("keypair");
		crypto_sign_keypair_low_stack(pk, sk);
		sp_end(&sc);
		// print_bytes_dec("pk", pk, CRYPTO_PUBLICKEYBYTES);

		// print_bytes_dec("SK", sk, CRYPTO_SECRETKEYBYTES);

		hal_send_str("KeyPair OK!\n");

		for (int i = 0; i < 1; i++)
		{

			sp_scope_t sp_sign = sp_begin("sign");
			crypto_sign_signature(sig, &siglen, m, mlen, sk);
			sp_end(&sp_sign);

			sp_scope_t sp_verify = sp_begin("verify");
			int res = crypto_sign_verify(sig, siglen, m, mlen, pk);
			sp_end(&sp_verify);

			// TODO: If verification fails, trap here (again: no I/O)
			if (res != 0)
			{
				hal_send_str("Signature Fail\n");
				return 1;
			}
			hal_send_str("Signature OK!\n");

			m[0]++;
		}

		hal_send_str("====== END ======");
	}

	return 0;
}
*/
