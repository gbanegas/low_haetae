// clang-format off
#ifndef HAETAE_POLY_H
#define HAETAE_POLY_H

#include "params.h"
#include "reduce.h"
#include "sampler.h"
#include <stdint.h>
#include <stddef.h>

typedef struct {
    int32_t coeffs[N];
} poly;

#if FROZEN_A >= 1
typedef struct {
    uint16_t coeffs[N];
} poly_frozen;

/* Multiply each coefficient by 2 modulo q (kept in uint16_t range). */
#define poly_frozen_double HAETAE_NAMESPACE(poly_frozen_double)
void poly_frozen_double(poly_frozen *a);

#define poly_uniform_frozen HAETAE_NAMESPACE(poly_uniform_frozen)
void poly_uniform_frozen(poly_frozen *a, const uint8_t seed[SEEDBYTES], uint16_t nonce);

#define poly_pointwise_montgomery_mixed_addto HAETAE_NAMESPACE(poly_pointwise_montgomery_mixed_addto)
void poly_pointwise_montgomery_mixed_addto(poly *c, const poly_frozen *a, const poly *b);

#define poly_pointwise_montgomery_mixed HAETAE_NAMESPACE(poly_pointwise_montgomery_mixed)
void poly_pointwise_montgomery_mixed(poly *c, const poly_frozen *a, const poly *b);
#endif

#if FROZEN_A >= 2
#define poly_uniform_frozen_pointwise_montgomery_set HAETAE_NAMESPACE(poly_uniform_frozen_pointwise_montgomery_set)
void poly_uniform_frozen_pointwise_montgomery_set(poly *dest, const uint8_t seed[SEEDBYTES], uint16_t nonce, const poly *src_hat);

#define poly_uniform_frozen_pointwise_montgomery_addto HAETAE_NAMESPACE(poly_uniform_frozen_pointwise_montgomery_addto)
void poly_uniform_frozen_pointwise_montgomery_addto(poly *acc, const uint8_t seed[SEEDBYTES], uint16_t nonce, const poly *src_hat);
#endif

#define poly_add HAETAE_NAMESPACE(poly_add)
void poly_add(poly *c, const poly *a, const poly *b);
#define poly_sub HAETAE_NAMESPACE(poly_sub)
void poly_sub(poly *c, const poly *a, const poly *b);
#define poly_pointwise_montgomery HAETAE_NAMESPACE(poly_pointwise_montgomery)
void poly_pointwise_montgomery(poly *c, const poly *a, const poly *b);

/*
 * In-place accumulate variant:
 *   c[i] += a[i] * b[i] * 2^{-32} (Montgomery reduce), for all i.
 *
 * This avoids allocating a temporary poly in *_pointwise_acc_* routines.
 */
#define poly_pointwise_montgomery_addto HAETAE_NAMESPACE(poly_pointwise_montgomery_addto)
void poly_pointwise_montgomery_addto(poly *c, const poly *a, const poly *b);

#define poly_reduce2q HAETAE_NAMESPACE(poly_reduce2q)
void poly_reduce2q(poly *a);
#define poly_freeze2q HAETAE_NAMESPACE(poly_freeze2q)
void poly_freeze2q(poly *a);
#define poly_freeze HAETAE_NAMESPACE(poly_freeze)
void poly_freeze(poly *a);

#define poly_highbits HAETAE_NAMESPACE(poly_highbits)
void poly_highbits(poly *a2, const poly *a);
#define poly_lowbits HAETAE_NAMESPACE(poly_lowbits)
void poly_lowbits(poly *a1, const poly *a);
#define poly_compose HAETAE_NAMESPACE(poly_compose)
void poly_compose(poly *a, const poly *ha, const poly *la);
#define poly_lsb HAETAE_NAMESPACE(poly_lsb)
void poly_lsb(poly *a0, const poly *a);

#define poly_uniform HAETAE_NAMESPACE(poly_uniform)
void poly_uniform(poly *a, const uint8_t seed[SEEDBYTES], uint16_t nonce);
#define poly_uniform_eta HAETAE_NAMESPACE(poly_uniform_eta)
void poly_uniform_eta(poly *a, const uint8_t seed[CRHBYTES], uint16_t nonce);
#define poly_uniform_add HAETAE_NAMESPACE(poly_uniform_add)
void poly_uniform_add(poly *a, const uint8_t seed[SEEDBYTES], uint16_t nonce);
#define poly_uniform_eta_add HAETAE_NAMESPACE(poly_uniform_eta_add)
void poly_uniform_eta_add(poly *a, const uint8_t seed[CRHBYTES], uint16_t nonce);
#define poly_challenge HAETAE_NAMESPACE(poly_challenge)
void poly_challenge(poly *c, const uint8_t highbits_lsb[POLYVECK_HIGHBITS_PACKEDBYTES + POLYC_PACKEDBYTES], const uint8_t mu[SEEDBYTES]);

#define poly_decomposed_pack HAETAE_NAMESPACE(poly_decomposed_pack)
void poly_decomposed_pack(uint8_t *buf, const poly *a);
#define poly_decomposed_unpack HAETAE_NAMESPACE(poly_decomposed_unpack)
void poly_decomposed_unpack(poly *a, const uint8_t *buf);

#define poly_pack_highbits HAETAE_NAMESPACE(poly_pack_highbits)
void poly_pack_highbits(uint8_t *buf, const poly *a);

#define poly_pack_lsb HAETAE_NAMESPACE(poly_pack_lsb)
void poly_pack_lsb(uint8_t *buf, const poly *a);

#define polyq_pack HAETAE_NAMESPACE(polyq_pack)
void polyq_pack(uint8_t *r, const poly *a);
#define polyq_unpack HAETAE_NAMESPACE(polyq_unpack)
void polyq_unpack(poly *r, const uint8_t *a);

#define polyeta_pack HAETAE_NAMESPACE(polyeta_pack)
void polyeta_pack(uint8_t *r, const poly *a);
#define polyeta_unpack HAETAE_NAMESPACE(polyeta_unpack)
void polyeta_unpack(poly *r, const uint8_t *a);
#define poly2eta_pack HAETAE_NAMESPACE(poly2eta_pack)
void poly2eta_pack(uint8_t *r, const poly *a);
#define poly2eta_unpack HAETAE_NAMESPACE(poly2eta_unpack)
void poly2eta_unpack(poly *r, const uint8_t *a);

#define poly_fromcrt HAETAE_NAMESPACE(poly_fromcrt)
void poly_fromcrt(poly *w, const poly *u, const poly *v);
#define poly_fromcrt0 HAETAE_NAMESPACE(poly_fromcrt0)
void poly_fromcrt0(poly *w, const poly *u);

#define poly_ntt HAETAE_NAMESPACE(poly_ntt)
void poly_ntt(poly *a);
#define poly_invntt_tomont HAETAE_NAMESPACE(poly_invntt_tomont)
void poly_invntt_tomont(poly *a);


#define poly_decompose_vk HAETAE_NAMESPACE(poly_decompose_vk)
void poly_decompose_vk(poly *v0, poly *v);

#define poly_decompose_stream_pack HAETAE_NAMESPACE(poly_decompose_stream_pack)
void poly_decompose_stream_pack(uint8_t pk[CRYPTO_PUBLICKEYBYTES], poly *v, size_t count);

#endif
// clang-format on
