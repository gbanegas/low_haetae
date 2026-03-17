#include "sign.h"
#include "packing.h"
#include "params.h"
#include "decompose.h"
#include "poly.h"
#include "polyfix.h"
#include "poly_sparsemul.h"
#include "polymat.h"
#include "polyvec.h"
#include "randombytes.h"
#include "symmetric.h"
#include "stack_profile.h"
#include <inttypes.h>
#include <stdint.h>
#include <string.h>


/* Optional: move large temporary vectors to BSS to fit small-stack profiles.
 * This is not thread-safe, but signing/verification on MCUs is typically single-threaded.
 */
#ifndef BRS_BSS_WORKSPACE
#define BRS_BSS_WORKSPACE 1
#endif

#if BRS_BSS_WORKSPACE
typedef union {
    polyfixvecl y1;
    polyfixvecl z1;
} brs_ws_v1_t;

typedef union {
    polyfixveck y2;
    polyfixveck z2;
} brs_ws_v2_t;

static struct {
    brs_ws_v1_t v1;
    brs_ws_v2_t v2;
    polyveck Ay;
} brs_sign_ws;
#endif



/*******************************************************************************
 * Name:        crypto_sign_signature
 *
 * Description: Computes signature.
 *
 * Arguments:   - uint8_t *sig:   pointer to output signature (of length
 *                                CRYPTO_BYTES)
 *              - size_t *siglen: pointer to output length of signature
 *              - uint8_t *m:     pointer to message to be signed
 *              - size_t mlen:    length of message
 *              - uint8_t *sk:    pointer to bit-packed secret key
 *
 * Returns 0 (success)
 **************************************************/

/* -------------------------------------------------------------------------
 * Low-stack signing strategy:
 *   - Pass A: sample y, compute c = H(HighBits(A*round(y)), LSB(round(y0)), mu)
 *   - Pass B: re-sample the same y (same counter_start), compute z = y + (-1)^b * LN * (c*s)
 *             using sparse-c multiplication (no NTT), and apply rejection tests.
 *   - Pass C: re-sample y again (same counter_start), recompute Ay = A*round(y) mod 2q,
 *             build the hint using (Ay, z2), and pack the signature.
 *
 * This avoids keeping Ay or the full secret vectors alive across the rejection loop.
 * ------------------------------------------------------------------------- */

#ifndef SIGN_RECOMPUTE_Y
#define SIGN_RECOMPUTE_Y 1
#endif

#if SIGN_RECOMPUTE_Y

/* Forward declarations (avoid implicit int / conflicting linkage on embedded builds). */
static __attribute__((noinline)) uint16_t
sample_y(polyfixvecl *y1, polyfixveck *y2, uint8_t *b,
         const uint8_t seedbuf[CRHBYTES], uint16_t counter_start);

static __attribute__((noinline)) void
compute_Ay_roundy_stream(polyveck *Ay,
                         const uint8_t *pkbytes,
                         polyvecl *z1rnd, polyveck *z2rnd,
                         const uint8_t z1rnd0_parity[N/8]);

/*
 * Ultra-low-stack recomputation path for the final hint/packing step.
 *
 * Goal: avoid keeping large temporary vectors alive in the same stack frame.
 * We split Pass C into three noinline helpers so that:
 *   - the (y1,y2,y1rnd,y2rnd) workspace exists only while rebuilding Ay;
 *   - the (z2rnd,htmp) workspace exists only while building the hint;
 *   - the (z1rnd) workspace exists only while packing the signature.
 *
 * This reduces the peak stack because these buffers no longer overlap.
 */

#if !BRS_BSS_WORKSPACE

static __attribute__((noinline)) void
recompute_Ay_from_seed(polyveck *Ay,
                       const uint8_t seedbuf[CRHBYTES],
                       uint16_t counter_start,
                       const uint8_t *pkbytes) {
    polyfixvecl y1;
    polyfixveck y2;
    uint8_t btmp = 0;
    (void)sample_y(&y1, &y2, &btmp, seedbuf, counter_start);

    /* Round in-place to avoid extra (y1rnd,y2rnd) temporaries. */
    polyfixvecl_round_inplace(&y1);
    polyfixveck_round_inplace(&y2);

    uint8_t y1rnd0_parity[N/8];
    memset(y1rnd0_parity, 0, sizeof y1rnd0_parity);

    for (unsigned i = 0; i < N; i++) {
        const uint8_t bit = (uint8_t)(y1.vec[0].coeffs[i] & 1);
        y1rnd0_parity[i >> 3] |= (uint8_t)(bit << (i & 7));
    }

    compute_Ay_roundy_stream(Ay, pkbytes,
                             &y1,
                             &y2,
                             y1rnd0_parity);
}

#endif

static __attribute__((noinline)) void
compute_hint_inplace(polyveck *Ay,
                     const polyfixveck *z2fix) {

    /* Overwrite Ay with h = HighBits(Ay) - HighBits(Ay - 2*round(z2)) (mod 2*ALPHA).
     *
     * BRS-style: avoid allocating (htmp,z2rnd) as full polyveck temporaries.
     * We do it row-by-row with a single poly scratch.
     */
    poly tmp;

    for (unsigned i = 0; i < K; i++) {
        /* tmp = Ay_i - 2*round(z2_i) */
        polyfix_round(&tmp, &z2fix->vec[i]);
        for (unsigned j = 0; j < N; j++) {
            tmp.coeffs[j] = Ay->vec[i].coeffs[j] - 2 * tmp.coeffs[j];
        }
        poly_freeze2q(&tmp);

        /* tmp = HighBits(tmp) */
        for (unsigned j = 0; j < N; j++) {
            decompose_hint(&tmp.coeffs[j], tmp.coeffs[j]);
        }

        /* Ay_i = HighBits(Ay_i) */
        for (unsigned j = 0; j < N; j++) {
            decompose_hint(&Ay->vec[i].coeffs[j], Ay->vec[i].coeffs[j]);
        }

        /* Ay_i = Ay_i - tmp */
        for (unsigned j = 0; j < N; j++) {
            Ay->vec[i].coeffs[j] -= tmp.coeffs[j];
        }
    }

    /* Center / wrap negative values. */
    polyveck_caddDQ2ALPHA(Ay);
}

static __attribute__((noinline)) int
pack_signature_from_z1fix(uint8_t *sig,
                          const poly *c,
                          polyfixvecl *z1fix,
                          const polyveck *h) {
    /* Round in place and pack using the same backing memory (polyfix == poly layout). */
    polyfixvecl_round_inplace(z1fix);
    return pack_sig_z1rnd(sig, c, (polyvecl *)(void *)z1fix, h);
}


static __attribute__((noinline)) void
stream_pk_ctx(uint8_t rhoprime[SEEDBYTES], polyveck *col0, const uint8_t *pkbytes) {
#if D == 1
    /* Low-stack pk streaming for D=1:
     *   - pk layout: [rhoprime || b1]
     *   - col0 (NTT) is derived as: col0[i] = NTT( 2*(a_i - 2*b1_i) mod Q )
     *     where a_i = poly_uniform(rhoprime, nonce_a(i)).
     *
     * We compute each col0[i] in place to avoid allocating full polyveck temporaries.
     */
    poly tmp;
    memcpy(rhoprime, pkbytes, SEEDBYTES);

    const uint8_t *bptr = pkbytes + SEEDBYTES;
    for (unsigned i = 0; i < K; ++i) {
        /* col0[i] <- a_i */
        poly_uniform(&col0->vec[i], rhoprime, (K << 8) + M + i);
        /* tmp <- b1_i */
        polyq_unpack(&tmp, bptr + i * POLYQ_PACKEDBYTES);

        /* col0[i] = 2*(a_i - 2*b1_i) mod Q */
        for (unsigned j = 0; j < N; ++j) {
            /* x in [-(2Q-2), Q-1] */
            int32_t x = col0->vec[i].coeffs[j] - (tmp.coeffs[j] << 1);
            /* bring into [0, Q-1] with two conditional subtractions */
            x += (int32_t)(2 * Q);
            if (x >= Q) x -= Q;
            if (x >= Q) x -= Q;

            x <<= 1;
            if (x >= Q) x -= Q;
            col0->vec[i].coeffs[j] = x;
        }
        poly_ntt(&col0->vec[i]);
    }
#elif D == 0
    /* for D==0, pk stores the NTT-domain column directly */
    unpack_pk(col0, rhoprime, pkbytes);
#else
#error "Not yet implemented."
#endif
}

static __attribute__((noinline)) uint16_t
sample_y(polyfixvecl *y1, polyfixveck *y2, uint8_t *b, const uint8_t seedbuf[CRHBYTES],
         uint16_t counter_start) {
    return polyfixveclk_sample_hyperball(y1, y2, b, seedbuf, counter_start);
}

static __attribute__((noinline)) void
compute_Ay_roundy_stream(polyveck *Ay,
                         const uint8_t *pkbytes,
                         polyvecl *z1rnd, polyveck *z2rnd,
                         const uint8_t z1rnd0_parity[N/8]) {
    uint8_t rhoprime[SEEDBYTES];
    polyveck col0;

    stream_pk_ctx(rhoprime, &col0, pkbytes);

    /* A*round(y) mod q = A1*round(y1) + 2*round(y2) mod q */
    polyvecl_ntt(z1rnd);
    polymatkl_pointwise_montgomery_stream(Ay, rhoprime, &col0, z1rnd);
    polyveck_invntt_tomont(Ay);
    polyveck_double(z2rnd);
    polyveck_add(Ay, Ay, z2rnd);

    /* Recover A*round(y) mod 2q. For the first component we only need
     * the parity of round(y0), not the full polynomial. */
    for (unsigned i = 0; i < N; ++i) {
        const int32_t xq  = Ay->vec[0].coeffs[i];
        const uint8_t bit = (uint8_t)((z1rnd0_parity[i >> 3] >> (i & 7)) & 1);
        Ay->vec[0].coeffs[i] = xq + (Q & -((xq ^ (int32_t)bit) & 1));
    }
    for (unsigned i = 1; i < K; ++i) {
        poly_fromcrt0(&Ay->vec[i], &Ay->vec[i]);
    }
    polyveck_freeze2q(Ay);
}

#if (HAETAE_MODE == 5)
/* Dense-binary challenge helpers (Mode 5)
 *
 * Mode 5 uses a dense binary challenge c with coefficients in {0,1}. The
 * sparse representation (idx/sign) used for modes 2/3 is not applicable.
 *
 * We therefore multiply c*s in the integer negacyclic ring Z[x]/(x^N+1)
 * (no modular wrap; values remain far below q), which keeps memory low.
 */

static __attribute__((noinline)) void
poly_mul_c_binary(poly *out, const poly *cbin, const poly *s) {
    memset(out->coeffs, 0, sizeof(out->coeffs));
    for (unsigned i = 0; i < N; ++i) {
        if (cbin->coeffs[i] == 0) continue;
        for (unsigned j = 0; j < N; ++j) {
            int k = (int)j - (int)i;
            if (k >= 0) out->coeffs[j] += s->coeffs[k];
            else        out->coeffs[j] -= s->coeffs[k + N];
        }
    }
}

static __attribute__((noinline)) uint64_t
polyfix_addmul_dense_binary_inplace_acc2(polyfix *z,
                                        const poly *s,
                                        const poly *cbin,
                                        int32_t k,
                                        uint64_t acc,
                                        const uint64_t bound) {
    poly cs;
    poly_mul_c_binary(&cs, cbin, s);
    for (unsigned j = 0; j < N; ++j) {
        const int32_t y = z->coeffs[j];
        const int32_t t = (int32_t)((int64_t)k * (int64_t)cs.coeffs[j]);
        const int32_t znew = y + t;
        const int32_t u = y + (t << 1);
        if (acc < bound) {
            uint64_t term2 = (uint64_t)((int64_t)u * (int64_t)u);
            if (term2 >= bound - acc) acc = bound;
            else acc += term2;
        }
        z->coeffs[j] = znew;
    }
    return acc;
}
#endif

#if BRS_BSS_WORKSPACE
static __attribute__((noinline)) void
recompute_Ay_from_seed_ws(polyveck *Ay,
                          polyfixvecl *y1,
                          polyfixveck *y2,
                          const uint8_t seedbuf[CRHBYTES],
                          uint16_t counter_start,
                          const uint8_t *pkbytes) {
    uint8_t b_tmp = 0;
    uint8_t y0_parity[N / 8];

    memset(y0_parity, 0, sizeof(y0_parity));

    /* Deterministically resample y into the provided buffers. */
    (void)sample_y(y1, y2, &b_tmp, seedbuf, counter_start);

    /* Round in place and compute (A1*y1 + A2*y2) in 2q. */
    polyfixvecl_round_inplace(y1);
    polyfixveck_round_inplace(y2);

    /* Pack parity bits of round(y0). */
    for (unsigned int i = 0; i < N; i++) {
        const uint8_t bit = (uint8_t)(y1->vec[0].coeffs[i] & 1);
        y0_parity[i / 8] |= (uint8_t)(bit << (i % 8));
    }

    compute_Ay_roundy_stream(Ay, pkbytes,
                             y1,
                             y2,
                             y0_parity);
}
#endif

static __attribute__((noinline)) void
passA_compute_challenge(poly *c, uint8_t *b_out,
                        uint16_t *counter_io,
                        const uint8_t seedbuf[CRHBYTES],
                        const uint8_t mu[CRHBYTES],
                        const uint8_t *pkbytes) {

    uint16_t counter_start = *counter_io;

    /* Challenge input = (HighBits(A*round(y) mod 2q), LSB(round(y0))).
     * BRS-style reduction: avoid keeping both fixed-point y and rounded y,
     * and avoid extra temporaries (lsb poly, highbits polyveck).
     */
    uint8_t buf[POLYVECK_HIGHBITS_PACKEDBYTES + POLYC_PACKEDBYTES] = {0};
#if BRS_BSS_WORKSPACE
    polyfixvecl *y1 = &brs_sign_ws.v1.y1;
    polyfixveck *y2 = &brs_sign_ws.v2.y2;
    polyveck    *Ay = &brs_sign_ws.Ay;
#else
    polyfixvecl y1_local;
    polyfixveck y2_local;
    polyveck    Ay_local;
    polyfixvecl *y1 = &y1_local;
    polyfixveck *y2 = &y2_local;
    polyveck    *Ay = &Ay_local;
#endif
    uint8_t y0_parity[N / 8];
    memset(y0_parity, 0, sizeof y0_parity);

    /* 1) sample y */
    *counter_io = sample_y(y1, y2, b_out, seedbuf, counter_start);

    /* 2) round y in-place (safe: same layout as poly/polyvec) */
    polyfixvecl_round_inplace(y1);
    polyfixveck_round_inplace(y2);

    /* 3) parity bits of round(y0) and LSB(round(y0)) packing */
    for (unsigned i = 0; i < N; i++) {
        const uint8_t bit = (uint8_t)(y1->vec[0].coeffs[i] & 1);
        y0_parity[i >> 3] |= (uint8_t)(bit << (i & 7));
        buf[POLYVECK_HIGHBITS_PACKEDBYTES + (i >> 3)] |= (uint8_t)(bit << (i & 7));
    }
    /* 4) Ay = A*round(y) mod 2q */
    compute_Ay_roundy_stream(Ay, pkbytes,
                             y1,
                             y2,
                             y0_parity);

    /* 5) HighBits(Ay) and challenge c = H(HighBits || LSB || mu) */
    polyveck_highbits_hint(Ay, Ay); /* in-place */
    polyveck_pack_highbits(buf, Ay);
    poly_challenge(c, buf, mu);
}

static __attribute__((noinline)) int
passB_compute_z_and_check(polyfixvecl *z1fix, polyfixveck *z2fix,
                          const poly *c, uint8_t b,
                          const uint8_t seedbuf[CRHBYTES],
                          uint16_t counter_start,
                          const uint8_t *sk_s1_packed,
                          const uint8_t *sk_s2_packed) {
    /* Re-sample the same y and overwrite it with z in-place */
    uint8_t b2 = 0;
    (void)b2;

    const uint64_t bound0 = (uint64_t)B0SQ * (uint64_t)LN * (uint64_t)LN;
    uint64_t acc2 = 0;
    const int32_t sgn = (b & 1) ? -1 : 1;

    /*
     * Challenge multiplication in Pass B
     * -------------------------------
     * Modes 2/3 use a sparse challenge (weight TAU), so we can multiply by c
     * using the sparse index/sign representation (BRS-style, no NTT).
     *
     * Mode 5 uses a dense binary challenge (0/1 coefficients). The sparse
     * representation is not applicable, so we switch to a low-stack dense
     * negacyclic product in the integer ring Z[x]/(x^N+1).
     */
#if (HAETAE_MODE == 2) || (HAETAE_MODE == 3)
    uint16_t idx[TAU];
    int8_t   sgns[TAU];
    unsigned w = poly_c_to_sparse(idx, sgns, TAU, c);
#endif

    /* 1) sample y into (z1fix,z2fix) */
    (void)sample_y(z1fix, z2fix, &b2, seedbuf, counter_start);

    /* 2) z1[0] = y1[0] + sgn*LN*c ; also accumulate acc2 for 2z - y */
    for (unsigned i = 0; i < N; ++i) {
        const int32_t y = z1fix->vec[0].coeffs[i];
        const int32_t t = (int32_t)((int64_t)sgn * (int64_t)LN * (int64_t)c->coeffs[i]);
        const int32_t z = y + t;
        const int32_t u = y + (t << 1);

        if (acc2 < bound0) {
            uint64_t term2 = (uint64_t)((int64_t)u * (int64_t)u);
            if (term2 >= bound0 - acc2) acc2 = bound0;
            else acc2 += term2;
        }
        z1fix->vec[0].coeffs[i] = z;
    }

    /* 3) remaining z1 components */
    for (unsigned i = 1; i < L; ++i) {
        poly s;
        polyeta_unpack(&s, sk_s1_packed + (i - 1) * POLYETA_PACKEDBYTES);
#if (HAETAE_MODE == 5)
        acc2 = polyfix_addmul_dense_binary_inplace_acc2(&z1fix->vec[i], &s, c,
                                                        (int32_t)((int64_t)sgn * (int64_t)LN),
                                                        acc2, bound0);
#else
        acc2 = polyfix_addmul_sparse_inplace_acc2(&z1fix->vec[i], &s,
                                                  idx, sgns, w,
                                                  (int32_t)((int64_t)sgn * (int64_t)LN),
                                                  acc2, bound0);
#endif
    }

    /* 4) z2 components */
    for (unsigned i = 0; i < K; ++i) {
        poly s;
#if D == 1
        poly2eta_unpack(&s, sk_s2_packed + i * POLY2ETA_PACKEDBYTES);
#elif D == 0
        polyeta_unpack(&s, sk_s2_packed + i * POLYETA_PACKEDBYTES);
#else
#error "Not yet implemented."
#endif
#if (HAETAE_MODE == 5)
        acc2 = polyfix_addmul_dense_binary_inplace_acc2(&z2fix->vec[i], &s, c,
                                                        (int32_t)((int64_t)sgn * (int64_t)LN),
                                                        acc2, bound0);
#else
        acc2 = polyfix_addmul_sparse_inplace_acc2(&z2fix->vec[i], &s,
                                                  idx, sgns, w,
                                                  (int32_t)((int64_t)sgn * (int64_t)LN),
                                                  acc2, bound0);
#endif
    }

    /* 5) rejection tests */
    {
        const uint64_t acc1 = polyfixveclk_sqnorm2(z1fix, z2fix);
        uint64_t reject1 = ((uint64_t)B1SQ * (uint64_t)LN * (uint64_t)LN - acc1) >> 63;
        reject1 &= 1;

        uint64_t reject2 = (acc2 - bound0) >> 63;
        reject2 &= 1;
        reject2 &= (b & 0x2) >> 1; /* b' bit */

        return (int)(reject1 | reject2);
    }
}

static __attribute__((noinline)) int
passB_check_only(const poly *c,
                 uint8_t b,
                 const uint8_t seedbuf[CRHBYTES],
                 uint16_t counter_start,
                 const uint8_t *sk_s1_packed,
                 const uint8_t *sk_s2_packed) {
#if BRS_BSS_WORKSPACE
    return passB_compute_z_and_check(&brs_sign_ws.v1.z1, &brs_sign_ws.v2.z2, c, b,
                                    seedbuf, counter_start,
                                    sk_s1_packed, sk_s2_packed);
#else
    polyfixvecl z1fix;
    polyfixveck z2fix;
    return passB_compute_z_and_check(&z1fix, &z2fix, c, b,
                                    seedbuf, counter_start,
                                    sk_s1_packed, sk_s2_packed);
#endif
}

static __attribute__((noinline)) int
passC_recompute_z_hint_and_pack(uint8_t *sig, size_t *siglen,
                                const poly *c,
                                uint8_t b,
                                const uint8_t seedbuf[CRHBYTES],
                                uint16_t counter_start,
                                const uint8_t *pkbytes,
                                const uint8_t *sk_s1_packed,
                                const uint8_t *sk_s2_packed) {

#if BRS_BSS_WORKSPACE
    polyveck *Ay = &brs_sign_ws.Ay;

    /* Rebuild Ay = A*round(y) mod 2q from seedbuf/counter_start using BSS buffers. */
    recompute_Ay_from_seed_ws(Ay, &brs_sign_ws.v1.y1, &brs_sign_ws.v2.y2,
                             seedbuf, counter_start, pkbytes);

    /* Recompute z and re-run checks (deterministic for given seedbuf/counter_start). */
    if (passB_compute_z_and_check(&brs_sign_ws.v1.z1, &brs_sign_ws.v2.z2, c, b,
                                 seedbuf, counter_start,
                                 sk_s1_packed, sk_s2_packed)) {
        return 1;
    }

    /* Overwrite Ay with hint. */
    compute_hint_inplace(Ay, &brs_sign_ws.v2.z2);

    /* Pack signature using round(z1) and hint. */
    if (pack_signature_from_z1fix(sig, c, &brs_sign_ws.v1.z1, Ay)) {
        return 1;
    }

    *siglen = CRYPTO_BYTES;
    return 0;

#else
    polyfixvecl z1fix;
    polyfixveck z2fix;

    /* Recompute z and re-run checks (deterministic for given seedbuf/counter_start). */
    if (passB_compute_z_and_check(&z1fix, &z2fix, c, b,
                                 seedbuf, counter_start,
                                 sk_s1_packed, sk_s2_packed)) {
        return 1;
    }

    /* Rebuild Ay = A*round(y) mod 2q from seedbuf/counter_start. */
    polyveck Ay;
    recompute_Ay_from_seed(&Ay, seedbuf, counter_start, pkbytes);

    /* Overwrite Ay with hint. */
    compute_hint_inplace(&Ay, &z2fix);

    /* Pack signature using round(z1) and hint. */
    if (pack_signature_from_z1fix(sig, c, &z1fix, &Ay)) {
        return 1;
    }

    *siglen = CRYPTO_BYTES;
    return 0;
#endif
}

#endif /* SIGN_RECOMPUTE_Y */

int crypto_sign_signature(uint8_t *sig, size_t *siglen, const uint8_t *m,
                          size_t mlen, const uint8_t *sk) {

   

    /* ------------------------------------------------------------
     * 0) Derive mu and per-signature seedbuf
     * ------------------------------------------------------------ */
    uint8_t seedbuf[CRHBYTES] = {0};
    uint8_t key[SEEDBYTES] = {0};
    uint8_t mu[CRHBYTES] = {0};
    xof256_state state;

    xof256_absorbe_twice(&state, sk, CRYPTO_PUBLICKEYBYTES, m, mlen);
    xof256_squeeze(mu, CRHBYTES, &state);

#if STREAM_MATRIX
    /* secret key layout: [pk || s1 || s2 || key] */
    const uint8_t *pkbytes = sk;
    const uint8_t *skp = sk + CRYPTO_PUBLICKEYBYTES;
    const uint8_t *sk_s1_packed = skp;
    skp += M * POLYETA_PACKEDBYTES;
#if D == 1
    const uint8_t *sk_s2_packed = skp;
    skp += K * POLY2ETA_PACKEDBYTES;
#elif D == 0
    const uint8_t *sk_s2_packed = skp;
    skp += K * POLYETA_PACKEDBYTES;
#else
#error "Not yet implemented."
#endif
    memcpy(key, skp, SEEDBYTES);
#else
    /* non-stream fallback: unpack full sk (kept for compatibility) */
    polyvecl A1[K];
    polyvecm s1;
    polyveck s2;
    unpack_sk(A1, &s1, &s2, key, sk);
    (void)A1; (void)s1; (void)s2;
    const uint8_t *pkbytes = sk;
    const uint8_t *sk_s1_packed = NULL;
    const uint8_t *sk_s2_packed = NULL;
#endif

    xof256_absorbe_twice(&state, key, SEEDBYTES, mu, CRHBYTES);
    xof256_squeeze(seedbuf, CRHBYTES, &state);

    /* ------------------------------------------------------------
     * Rejection sampling loop
     * ------------------------------------------------------------ */
    uint16_t counter = 0;

    for (;;) {
        uint16_t counter_start = counter;
        uint8_t b = 0;
        poly c;

#if SIGN_RECOMPUTE_Y && STREAM_MATRIX
        
        passA_compute_challenge(&c, &b, &counter, seedbuf, mu, pkbytes);


        
        if (passB_check_only(&c, b, seedbuf, counter_start,
                             sk_s1_packed, sk_s2_packed)) {
    
            continue; /* reject */
        }

        if (passC_recompute_z_hint_and_pack(sig, siglen, &c, b,
                                            seedbuf, counter_start,
                                            pkbytes, sk_s1_packed, sk_s2_packed)) {
        
            continue; /* encoding failed => reject */
        }
     
        return 0;
#else
        /* If SIGN_RECOMPUTE_Y is disabled or STREAM_MATRIX is off, fall back to
         * the reference implementation in this file (not provided here). */
        (void)counter_start; (void)b; (void)c;
        return -1;
#endif
    }
}


/********************************************************************
 * Name:        crypto_sign
 *
 * Description: Compute signed message.
 *
 * Arguments:   - uint8_t *sm: pointer to output signed message (allocated
 *                             array with CRYPTO_BYTES + mlen bytes),
 *                             can be equal to m
 *              - size_t *smlen: pointer to output length of signed
 *                               message
 *              - const uint8_t *m: pointer to message to be signed
 *              - size_t mlen: length of message
 *              - const uint8_t *sk: pointer to bit-packed secret key
 *
 * Returns 0 (success)
 **************************************************/
int crypto_sign_sign(uint8_t *sm, size_t *smlen, const uint8_t *m, size_t mlen,
                     const uint8_t *sk) {
    size_t i;

    for (i = 0; i < mlen; ++i)
        sm[CRYPTO_BYTES + mlen - 1 - i] = m[mlen - 1 - i];
    crypto_sign_signature(sm, smlen, sm + CRYPTO_BYTES, mlen, sk);
    *smlen += mlen;
    return 0;
}
