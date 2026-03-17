#include <inttypes.h>
#include <stdint.h>
#include <string.h>

#include "config.h"
#include "packing.h"
#include "params.h"
#include "poly.h"
#include "polyfix.h"
#include "polymat.h"
#include "polyvec.h"
#include "randombytes.h"
#include "stack_profile.h"
#include "symmetric.h"
#include "utils.h"

#define crypto_sign_keypair HAETAE_NAMESPACE(keypair)
int crypto_sign_keypair(uint8_t *pk, uint8_t *sk);

#define crypto_sign_keypair_low_stack HAETAE_NAMESPACE(keypair_low)
int crypto_sign_keypair_low_stack(uint8_t *pk, uint8_t *sk);

#define crypto_sign_keypair_streamed_v1 HAETAE_NAMESPACE(crypto_sign_keypair_streamed_v1)
int crypto_sign_keypair_streamed_v1(uint8_t *pk, uint8_t *sk);

#define crypto_sign_keypair_streamed_v2 HAETAE_NAMESPACE(crypto_sign_keypair_streamed_v2)
int crypto_sign_keypair_streamed_v2(uint8_t *pk, uint8_t *sk);
