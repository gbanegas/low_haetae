#ifndef HAETAE_SIGN_H
#define HAETAE_SIGN_H

#include <inttypes.h>
#include <stdint.h>
#include <string.h>
#include "packing.h"
#include "params.h"
#include "poly.h"
#include "polyfix.h"
#include "polymat.h"
#include "polyvec.h"
#include "randombytes.h"
#include "symmetric.h"
#include "config.h"

#define crypto_sign_signature HAETAE_NAMESPACE(signature)
int crypto_sign_signature(uint8_t *sig, size_t *siglen, const uint8_t *m,
                          size_t mlen, const uint8_t *sk);

#define crypto_sign_sign HAETAE_NAMESPACE(sign)
int crypto_sign_sign(uint8_t *sm, size_t *smlen, const uint8_t *m, size_t mlen,
                     const uint8_t *sk);


#endif // HAETAE_SIGN_H