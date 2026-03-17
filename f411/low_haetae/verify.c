#include "verify.h"

/* Optional debug helpers for investigating signature failures.
 * Enable with -DDEBUG_MODE5_TRACE=1.
 */
#if defined(DEBUG_MODE5_TRACE) && DEBUG_MODE5_TRACE
#include <stdio.h>

static void dbg_hex(const char *label, const uint8_t *x, size_t n)
{
    printf("[DBG] %s (%zu): ", label, n);
    for (size_t i = 0; i < n; i++)
        printf("%02x", x[i]);
    printf("\n");
}

static void dbg_pack_c(uint8_t out[POLYC_PACKEDBYTES], const poly *c)
{
    for (size_t i = 0; i < POLYC_PACKEDBYTES; i++) out[i] = 0;
    for (size_t i = 0; i < N; i++) {
        const uint8_t bit = (uint8_t)(c->coeffs[i] & 1);
        out[i >> 3] |= (uint8_t)(bit << (i & 7));
    }
}
#endif

/* -------------------------- helpers (only used for STREAM_MATRIX && (D==0||D==1)) -------------------------- */
#if STREAM_MATRIX && ((D == 1) || (D == 0))
static inline uint8_t get_bit(const uint8_t *a, size_t i) {
    return (uint8_t)((a[i >> 3] >> (i & 7)) & 1u);
}

static inline void set_bit(uint8_t *a, size_t i, uint8_t b) {
    a[i >> 3] |= (uint8_t)((b & 1u) << (i & 7));
}

/*
 * poly_fromcrt() variant where the CRT "mod 2" component is provided
 * as packed bits (N/8 bytes).
 */
static void poly_fromcrt_bits(poly *w, const poly *u,
                              const uint8_t vbits[POLYC_PACKEDBYTES]) {
    for (size_t i = 0; i < N; i++) {
        int32_t xq = u->coeffs[i];
        int32_t x2 = (int32_t)get_bit(vbits, i);
        w->coeffs[i] = xq + (Q & -((xq ^ x2) & 1));
    }
}

/* Add a packed 0/1 polynomial (bits) into a poly. */
static void poly_add_bits(poly *a, const uint8_t bits[POLYC_PACKEDBYTES]) {
    for (size_t i = 0; i < N; i++) {
        a->coeffs[i] += (int32_t)get_bit(bits, i);
    }
}

/* Per-polynomial helpers to avoid polyveck temporaries in low-SRAM paths. */
static void poly_highbits_hint_one(poly *w, const poly *v) {
    for (size_t j = 0; j < N; j++) {
        decompose_hint(&w->coeffs[j], v->coeffs[j]);
    }
}

static void poly_csubDQ2ALPHA_one(poly *v) {
    for (size_t j = 0; j < N; j++) {
        v->coeffs[j] -= ~((v->coeffs[j] - (DQ - 2) / ALPHA_HINT) >> 31) &
                        ((DQ - 2) / ALPHA_HINT);
    }
}

static void poly_div2_one(poly *v) {
    for (size_t j = 0; j < N; j++) {
        v->coeffs[j] >>= 1;
    }
}

static uint64_t poly_sqnorm2_one(const poly *b) {
    uint64_t ret = 0;
    for (size_t j = 0; j < N; j++) {
        ret += (uint64_t)b->coeffs[j] * (uint64_t)b->coeffs[j];
    }
    return ret;
}
#endif

/*************************************************
 * Name:        crypto_sign_verify
 *
 * Description: Verifies signature.
 *
 * Arguments:   - const uint8_t *sig: pointer to signature
 *              - size_t siglen: length of signature
 *              - const uint8_t *m: pointer to message
 *              - size_t mlen: length of message
 *              - const uint8_t *pk: pointer to bit-packed public key
 *
 * Returns 0 if signature could be verified correctly and -1 otherwise
 **************************************************/
int crypto_sign_verify(const uint8_t *sig, size_t siglen, const uint8_t *m,
                       size_t mlen, const uint8_t *pk) {

    if (siglen != CRYPTO_BYTES) {
        return -1;
    }

#if STREAM_MATRIX && ((D == 1) || (D == 0))
    /*
     * Low-SRAM verification (BRS-style) for D==0/1:
     *   - decode HB(z1) to int8_t
     *   - compose/NTT z1 column-by-column
     *   - keep w' as packed bits
     *   - stream A1 multiplication
     */
    uint8_t buf[POLYVECK_HIGHBITS_PACKEDBYTES + POLYC_PACKEDBYTES] = {0};
    uint8_t rhoprime[SEEDBYTES] = {0}, mu[SEEDBYTES];
    int8_t hb_z1_i8[N * L];
    uint8_t wprime_bits[POLYC_PACKEDBYTES] = {0};

#if D == 0
    /* D==0 path: keep only per-row temporaries. */
    uint16_t h_u16[N * K];
    poly col0; /* pk row polynomial: NTT(-2b) */
    poly trow; /* running accumulator for one output row */
    poly z;    /* reused: z1 column / cprime */
#else
    /* D==1 path: keep the previous polyveck-based implementation. */
    polyveck col; /* b -> col0 (NTT) -> decoded h */
    polyveck t;
    polyveck w2;
    poly z; /* reused: z1 column / cprime */
#endif

#if FROZEN_A >= 1
    poly_frozen af;
#else
    poly a;
#endif

    uint64_t sqnorm2 = 0;
    xof256_state state;

    /* ---- parse signature (robust) ---- */
    const uint8_t *sp = sig;
    const uint8_t *c_bytes = sp;
    sp += POLYC_PACKEDBYTES;

    const uint8_t *lowbits_ptr = sp;
    sp += (size_t)L * (size_t)N;

    uint16_t size_enc_hb_z1 = (uint16_t)sp[0] + (uint16_t)BASE_ENC_HB_Z1;
    uint16_t size_enc_h = (uint16_t)sp[1] + (uint16_t)BASE_ENC_H;
    sp += 2;

    const size_t used = (size_t)POLYC_PACKEDBYTES + (size_t)L * (size_t)N + 2u +
                        (size_t)size_enc_hb_z1 + (size_t)size_enc_h;
    if (used > CRYPTO_BYTES) {
        return -1;
    }

    const uint8_t *hb_ptr = sp;
    const uint8_t *h_ptr = sp + size_enc_hb_z1;
    const uint8_t *pad_ptr = h_ptr + size_enc_h;

    if (decode_hb_z1_i8(hb_z1_i8, hb_ptr, size_enc_hb_z1)) {
        return -1;
    }

    /* Zero padding check. */
    for (size_t i = 0; i < (CRYPTO_BYTES - used); i++) {
        if (pad_ptr[i] != 0) {
            return -1;
        }
    }

    /* === D==0: row-streamed verification to minimize stack === */
#if D == 0
    /* Decode h directly to uint16_t. */
    if (decode_h_u16(h_u16, h_ptr, size_enc_h)) {
        return -1;
    }

    /* Extract seed from pk (first SEEDBYTES bytes). */
    memcpy(rhoprime, pk, SEEDBYTES);

    /*
     * Pre-pass: compute w' bits (only depends on ell==0) and ||z1||^2.
     * We will re-compose columns later for the (slow) row-wise NTT products.
     */
    for (size_t ell = 0; ell < L; ell++) {
        const size_t off = ell * (size_t)N;
        for (size_t j = 0; j < N; j++) {
            const int32_t hb = (int32_t)hb_z1_i8[off + j];
            const int32_t lo = (int32_t)(int8_t)lowbits_ptr[off + j];
            const int32_t zc = (hb << 8) + lo;
            sqnorm2 += (uint64_t)((int64_t)zc * (int64_t)zc);

            if (ell == 0) {
                const uint8_t cbit = get_bit(c_bytes, j);
                const uint8_t wbit = (uint8_t)((zc - (int32_t)cbit) & 1);
                set_bit(wprime_bits, j, wbit);
            }
        }
    }

    /*
     * Row-wise streamed multiplication: recompute each output row trow of
     * t = A1 * round(z1) without keeping polyveck temporaries.
     *
     * This trades speed for SRAM: each row re-composes and NTTs all L columns.
     */
    uint64_t z2_sqnorm2 = 0;
    for (size_t r = 0; r < K; r++) {
        /* Load pk row polynomial (NTT(-2b)) for ell==0. */
        polyq_unpack(&col0, pk + SEEDBYTES + r * (size_t)POLYQ_PACKEDBYTES);

        /* Clear accumulator (NTT domain). */
        memset(&trow, 0, sizeof(poly));

        for (size_t ell = 0; ell < L; ell++) {
            const size_t off = ell * (size_t)N;
            for (size_t j = 0; j < N; j++) {
                const int32_t hb = (int32_t)hb_z1_i8[off + j];
                const int32_t lo = (int32_t)(int8_t)lowbits_ptr[off + j];
                z.coeffs[j] = (hb << 8) + lo;
            }

            poly_ntt(&z);

            if (ell == 0) {
                poly_pointwise_montgomery(&trow, &col0, &z);
            } else {
                const uint16_t nonce = (uint16_t)((r << 8) + (uint16_t)(ell - 1));
#if FROZEN_A >= 1
                poly_uniform_frozen(&af, rhoprime, nonce);
                poly_frozen_double(&af);
                poly_pointwise_montgomery_mixed_addto(&trow, &af, &z);
#else
                poly_uniform(&a, rhoprime, nonce);
                poly_double(&a);
                poly_pointwise_montgomery_addto(&trow, &a, &z);
#endif
            }
        }

        /* Back to standard domain. */
        poly_invntt_tomont(&trow);

        /* Recover mod 2q using CRT (first row uses w'). */
        if (r == 0) {
            poly_fromcrt_bits(&trow, &trow, wprime_bits);
        } else {
            poly_fromcrt0(&trow, &trow);
        }
        poly_freeze2q(&trow);

        /* w_row = HighBits(t_row) + h_row. Use z as workspace. */
        poly_highbits_hint_one(&z, &trow);
        for (size_t j = 0; j < N; j++) {
            z.coeffs[j] += (int32_t)h_u16[r * (size_t)N + j];
        }
        poly_csubDQ2ALPHA_one(&z);

        /* Pack transcript HighBits(w_row). */
        poly_pack_highbits(buf + r * (size_t)POLY_HIGHBITS_PACKEDBYTES, &z);

        /* z2_row = (alpha*w_row - t_row + w'_row)/2 mod q and accumulate ||z2||^2. */
        for (size_t j = 0; j < N; j++) {
            trow.coeffs[j] = (int32_t)(ALPHA_HINT * (int64_t)z.coeffs[j]) - trow.coeffs[j];
        }
        if (r == 0) {
            poly_add_bits(&trow, wprime_bits);
        }
        poly_reduce2q(&trow);
        poly_div2_one(&trow);

        z2_sqnorm2 += poly_sqnorm2_one(&trow);
        if (sqnorm2 + z2_sqnorm2 > B2SQ) {
            return -1;
        }
    }

    /* Append w' bits to the transcript and recompute c'. */
    memcpy(buf + POLYVECK_HIGHBITS_PACKEDBYTES, wprime_bits, POLYC_PACKEDBYTES);

    xof256_absorbe_twice(&state, pk, CRYPTO_PUBLICKEYBYTES, m, mlen);
    xof256_squeeze(mu, SEEDBYTES, &state);

    poly_challenge(&z, buf, mu);

    for (size_t j = 0; j < N; j++) {
        if (get_bit(c_bytes, j) != (uint8_t)z.coeffs[j]) {
            return -1;
        }
    }

    return 0;

#elif D == 1
    /* === D==1: previous implementation === */
    /* Unpack public key: b-vector and seed. */
    unpack_pk(&col, rhoprime, pk);

#if D == 1
    /* ---- build col0 = 2*(a - 2*b) in-place in 'col', then NTT(col0) ---- */
    for (size_t r = 0; r < K; r++) {
        const uint16_t nonce_a = (uint16_t)((K << 8) + M + (uint16_t)r);
#if FROZEN_A >= 1
        poly_uniform_frozen(&af, rhoprime, nonce_a);
        for (size_t j = 0; j < N; j++) {
            int32_t ai = (int32_t)af.coeffs[j];
            int32_t bi = col.vec[r].coeffs[j];
            col.vec[r].coeffs[j] = (ai << 1) - (bi << 2); /* 2a - 4b */
        }
#else
        poly_uniform(&a, rhoprime, nonce_a);
        for (size_t j = 0; j < N; j++) {
            int32_t ai = a.coeffs[j];
            int32_t bi = col.vec[r].coeffs[j];
            col.vec[r].coeffs[j] = (ai << 1) - (bi << 2); /* 2a - 4b */
        }
#endif
    }
    polyveck_ntt(&col);
#endif

    /* ---- streamed multiplication t = A1 * round(z1) in NTT domain ---- */
    for (size_t ell = 0; ell < L; ell++) {
        /* Compose this column of z1: z = HB*256 + LB. */
        const size_t off = ell * (size_t)N;
        for (size_t j = 0; j < N; j++) {
            const int32_t hb = (int32_t)hb_z1_i8[off + j];
            const int32_t lo = (int32_t)(int8_t)lowbits_ptr[off + j];
            const int32_t zc = (hb << 8) + lo;
            z.coeffs[j] = zc;
            sqnorm2 += (uint64_t)((int64_t)zc * (int64_t)zc);

            if (ell == 0) {
                const uint8_t cbit = get_bit(c_bytes, j);
                const uint8_t wbit = (uint8_t)((zc - (int32_t)cbit) & 1);
                set_bit(wprime_bits, j, wbit);
            }
        }

        poly_ntt(&z);

        if (ell == 0) {
            /* First column: use precomputed col0. */
            for (size_t r = 0; r < K; r++) {
                poly_pointwise_montgomery(&t.vec[r], &col.vec[r], &z);
            }
        } else {
            /* Other columns: 2*A0 entries streamed from seed. */
            for (size_t r = 0; r < K; r++) {
                const uint16_t nonce = (uint16_t)((r << 8) + (uint16_t)(ell - 1));
#if FROZEN_A >= 1
                poly_uniform_frozen(&af, rhoprime, nonce);
                poly_frozen_double(&af);
                poly_pointwise_montgomery_mixed_addto(&t.vec[r], &af, &z);
#else
                poly_uniform(&a, rhoprime, nonce);
                poly_double(&a);
                poly_pointwise_montgomery_addto(&t.vec[r], &a, &z);
#endif
            }
        }
    }

    /* Back to standard domain. */
    polyveck_invntt_tomont(&t);

    /* Recover mod 2q using the CRT (first row uses w'). */
    poly_fromcrt_bits(&t.vec[0], &t.vec[0], wprime_bits);
    for (size_t r = 1; r < K; r++) {
        poly_fromcrt0(&t.vec[r], &t.vec[r]);
    }
    polyveck_freeze2q(&t);

    /* Decode h into col (col no longer needed as col0). */
    if (decode_h(&col.vec[0].coeffs[0], h_ptr, size_enc_h)) {
        return -1;
    }

    /* Recover w = HighBits(t) + h. */
    polyveck_highbits_hint(&w2, &t);
    polyveck_add(&w2, &w2, &col);
    polyveck_csubDQ2ALPHA(&w2);

    /* Pack challenge transcript: HighBits(w) and w'. */
    polyveck_pack_highbits(buf, &w2);
    memcpy(buf + POLYVECK_HIGHBITS_PACKEDBYTES, wprime_bits, POLYC_PACKEDBYTES);

    /* Recover z2 = (alpha*w - t + w')/2 mod q, and check norm. */
    polyveck_mul_alpha(&w2, &w2);
    polyveck_sub(&w2, &w2, &t);
    poly_add_bits(&w2.vec[0], wprime_bits);
    polyveck_reduce2q(&w2);
    polyveck_div2(&w2);

    if (sqnorm2 + polyveck_sqnorm2(&w2) > B2SQ) {
        return -1;
    }

    /* Recompute c' and compare with packed c in the signature. */
    xof256_absorbe_twice(&state, pk, CRYPTO_PUBLICKEYBYTES, m, mlen);
    xof256_squeeze(mu, SEEDBYTES, &state);

    poly_challenge(&z, buf, mu);

    for (size_t j = 0; j < N; j++) {
        if (get_bit(c_bytes, j) != (uint8_t)z.coeffs[j]) {
            return -1;
        }
    }

    return 0;

#endif /* D */

#else
    /* Fallback reference path (kept for completeness). */
    unsigned int i;
    uint8_t buf[POLYVECK_HIGHBITS_PACKEDBYTES + POLYC_PACKEDBYTES] = {0};
    uint8_t rhoprime[SEEDBYTES] = {0}, mu[SEEDBYTES];
    uint64_t sqnorm2;
#if !STREAM_MATRIX
    polyvecl A1[K];
#endif
    polyvecl lowbits_z1, z1;
    polyveck b, highbits, h, z2, w;
#if D > 0
    polyveck a;
#endif
    poly c, cprime, wprime;

    xof256_state state;

    /* Unpack public key */
    unpack_pk(&b, rhoprime, pk);

    /* Unpack signature */
    if (unpack_sig(&c, &lowbits_z1, &z1, &h, sig)) {
        return -1;
    }

    /* Compose z1 out of HighBits(z1) and LowBits(z1) */
    for (i = 0; i < L; ++i) {
        poly_compose(&z1.vec[i], &z1.vec[i], &lowbits_z1.vec[i]);
    }

    /* Recover A1 */
#if !STREAM_MATRIX
    polymatkl_expand(A1, rhoprime);
    polymatkl_double(A1);
#if D == 1
    polyveck_expand(&a, rhoprime);
    polyveck_double(&b);
    polyveck_sub(&b, &a, &b);
    polyveck_double(&b);
    polyveck_ntt(&b);
#elif D == 0
    /* pk already stores the NTT-domain column */
#else
#error "Not yet implemented."
#endif
    for (i = 0; i < K; ++i) {
        A1[i].vec[0] = b.vec[i];
    }
#endif

    /* Compute z2 */
    sqnorm2 = polyvecl_sqnorm2(&z1);
    poly_sub(&wprime, &z1.vec[0], &c);
    poly_lsb(&wprime, &wprime);

    polyvecl_ntt(&z1);
#if STREAM_MATRIX
    {
        polyveck b1 = b;
#if D > 0
        polyveck_expand(&a, rhoprime);
#endif
#if D == 1
        polyveck_double(&b1);
        polyveck_sub(&b1, &a, &b1);
        polyveck_double(&b1);
        polyveck_ntt(&b1);
#elif D == 0
        /* pk already stores the NTT-domain column */
#else
#error "Not yet implemented."
#endif
        polymatkl_pointwise_montgomery_stream(&highbits, rhoprime, &b1, &z1);
    }
#else
    polymatkl_pointwise_montgomery(&highbits, A1, &z1);
#endif
    polyveck_invntt_tomont(&highbits);

    polyveck_poly_fromcrt(&highbits, &highbits, &wprime);
    polyveck_freeze2q(&highbits);

    polyveck_highbits_hint(&w, &highbits);
    polyveck_add(&w, &w, &h);
    polyveck_csubDQ2ALPHA(&w);

    polyveck_mul_alpha(&z2, &w);
    polyveck_sub(&z2, &z2, &highbits);
    poly_add(&z2.vec[0], &z2.vec[0], &wprime);
    polyveck_reduce2q(&z2);
    polyveck_div2(&z2);

    if (sqnorm2 + polyveck_sqnorm2(&z2) > B2SQ) {
        return -1;
    }

    polyveck_pack_highbits(buf, &w);
    poly_pack_lsb(buf + POLYVECK_HIGHBITS_PACKEDBYTES, &wprime);

    xof256_absorbe_twice(&state, pk, CRYPTO_PUBLICKEYBYTES, m, mlen);
    xof256_squeeze(mu, SEEDBYTES, &state);

    poly_challenge(&cprime, buf, mu);

    for (i = 0; i < N; ++i) {
        if (c.coeffs[i] != cprime.coeffs[i]) {
#if defined(DEBUG_MODE5_TRACE) && DEBUG_MODE5_TRACE
            uint8_t c_bytes[POLYC_PACKEDBYTES];
            uint8_t cp_bytes[POLYC_PACKEDBYTES];
            dbg_pack_c(c_bytes, &c);
            dbg_pack_c(cp_bytes, &cprime);
            dbg_hex("mu", mu, 32);
            dbg_hex("c", c_bytes, POLYC_PACKEDBYTES);
            dbg_hex("c'", cp_bytes, POLYC_PACKEDBYTES);
            dbg_hex("buf[0..31]", buf, 32);
            dbg_hex("buf[last 32]", buf + POLYVECK_HIGHBITS_PACKEDBYTES, 32);
#endif
            return -1;
        }
    }
    return 0;
#endif
}

/*************************************************
 * Name:        crypto_sign_open
 *
 * Description: Verify signed message.
 *
 * Returns 0 if signed message could be verified correctly and -1 otherwise
 **************************************************/
int crypto_sign_open(uint8_t *m, size_t *mlen, const uint8_t *sm, size_t smlen,
                     const uint8_t *pk) {
    size_t i;

    if (smlen < CRYPTO_BYTES)
        goto badsig;

    *mlen = smlen - CRYPTO_BYTES;
    if (crypto_sign_verify(sm, CRYPTO_BYTES, sm + CRYPTO_BYTES, *mlen, pk))
        goto badsig;
    else {
        /* All good, copy msg, return 0 */
        for (i = 0; i < *mlen; ++i)
            m[i] = sm[CRYPTO_BYTES + i];
        return 0;
    }

badsig:
    /* Signature verification failed */
    *mlen = (size_t)-1;
    for (i = 0; i < smlen; ++i)
        m[i] = 0;

    return -1;
}
