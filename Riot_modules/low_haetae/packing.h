// clang-format off
#ifndef HAETAE_PACKING_H
#define HAETAE_PACKING_H

#include "params.h"
#include "polyvec.h"
#include <stdint.h>
#include <stddef.h>

#define pack_pk HAETAE_NAMESPACE(pack_pk)
void pack_pk(uint8_t pk[CRYPTO_PUBLICKEYBYTES], polyveck *b, const uint8_t seed[SEEDBYTES]);

/* Streamed packing helpers for low-stack keygen (D>0 path). */
#define pack_pk_stream HAETAE_NAMESPACE(pack_pk_stream)
void pack_pk_stream(uint8_t pk[CRYPTO_PUBLICKEYBYTES], poly *b, size_t count);

#define pack_pk_stream_remain HAETAE_NAMESPACE(pack_pk_stream_remain)
void pack_pk_stream_remain(uint8_t pk[CRYPTO_PUBLICKEYBYTES], const uint8_t seed[SEEDBYTES]);

#define unpack_pk HAETAE_NAMESPACE(unpack_pk)
void unpack_pk(polyveck *b, uint8_t seed[SEEDBYTES], const uint8_t pk[CRYPTO_PUBLICKEYBYTES]);

#define pack_sk HAETAE_NAMESPACE(pack_sk)
void pack_sk(uint8_t sk[CRYPTO_SECRETKEYBYTES], const uint8_t pk[CRYPTO_PUBLICKEYBYTES], const polyvecm *s0, const polyveck *s1, const uint8_t key[SEEDBYTES]);

/* Streamed packing helpers for low-stack keygen (D>0 path). */
#define pack_sk_stream_s1 HAETAE_NAMESPACE(pack_sk_stream_s1)
void pack_sk_stream_s1(uint8_t sk[CRYPTO_SECRETKEYBYTES], const poly s1, const size_t count);

#define pack_sk_stream_s2 HAETAE_NAMESPACE(pack_sk_stream_s2)
void pack_sk_stream_s2(uint8_t sk[CRYPTO_SECRETKEYBYTES], const poly s2, const size_t count);

#define pack_sk_stream_remain HAETAE_NAMESPACE(pack_sk_stream_remain)
void pack_sk_stream_remain(uint8_t sk[CRYPTO_SECRETKEYBYTES],
                           const uint8_t pk[CRYPTO_PUBLICKEYBYTES],
                           const uint8_t key[SEEDBYTES]);

#define unpack_sk HAETAE_NAMESPACE(unpack_sk)
void unpack_sk(polyvecl A[K], polyvecm *s0, polyveck *s1, uint8_t *key, const uint8_t sk[CRYPTO_SECRETKEYBYTES]);

#define pack_sig HAETAE_NAMESPACE(pack_sig)
int pack_sig(uint8_t sig[CRYPTO_BYTES], const poly *c, const polyvecl *lowbits_z1, const polyvecl *highbits_z1, const polyveck *h);

/*
 * Low-stack signing helper.
 *
 * Packs a signature sig = (c, LB(z1), len(x), len(y), x = Enc(HB(z1)), y = Enc(h))
 * while avoiding the temporary polyvecl buffers for LB(z1) and HB(z1).
 *
 * Input:  z1rnd is the rounded vector round(z1).
 * Effect: this function packs LB(round(z1)) directly to sig and then overwrites
 *         z1rnd in-place with HB(round(z1)) to run the encoder.
 */
#define pack_sig_z1rnd HAETAE_NAMESPACE(pack_sig_z1rnd)
int pack_sig_z1rnd(uint8_t sig[CRYPTO_BYTES], const poly *c, polyvecl *z1rnd,
                   const polyveck *h);

#define unpack_sig HAETAE_NAMESPACE(unpack_sig)
int unpack_sig(poly *c, polyvecl *lowbits_z1, polyvecl *highbits_z1, polyveck *h, const uint8_t sig[CRYPTO_BYTES]);

#endif //HAETAE_PACKING_H
// clang-format on
