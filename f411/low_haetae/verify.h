#include <inttypes.h>
#include <stdint.h>
#include <string.h>
#include "packing.h"
#include "params.h"
#include "poly.h"
#include "polyfix.h"
#include "polymat.h"
#include "polyvec.h"
#include "decompose.h"
#include "encoding.h"
#include "randombytes.h"
#include "symmetric.h"
#include "config.h"


#define crypto_sign_verify HAETAE_NAMESPACE(verify)
int crypto_sign_verify(const uint8_t *sig, size_t siglen, const uint8_t *m,
                       size_t mlen, const uint8_t *pk);


#define crypto_sign_open HAETAE_NAMESPACE(open)
int crypto_sign_open(uint8_t *m, size_t *mlen, const uint8_t *sm, size_t smlen,
                     const uint8_t *pk);

