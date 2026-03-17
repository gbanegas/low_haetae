#include "polymat.h"
#include "params.h"
#include "poly.h"
#include "polyvec.h"
#include <stdint.h>

/*************************************************
 * Name:        polymat_expand
 *
 * Description: Implementation of ExpandA. Generates matrix A with uniformly
 *              random coefficients a_{i,j} by performing rejection
 *              sampling on the output stream of SHAKE128(rho|j|i)
 *              or AES256CTR(rho,j|i).
 *
 * Arguments:   - polyvecm mat[K]: output matrix k \times m
 *              - const uint8_t rho[]: byte array containing seed rho
 **************************************************/
void polymatkl_expand(polyvecl mat[K], const uint8_t rho[SEEDBYTES]) {
    unsigned int i, j;

    for (i = 0; i < K; ++i)
        for (j = 0; j < M; ++j)
            poly_uniform(&mat[i].vec[j + 1], rho, (i << 8) + j);
}

/*************************************************
 * Name:        polymat_expand
 *
 * Description: Implementation of ExpandA. Generates matrix A with uniformly
 *              random coefficients a_{i,j} by performing rejection
 *              sampling on the output stream of SHAKE128(rho|j|i)
 *              or AES256CTR(rho,j|i).
 *
 * Arguments:   - polyvecm mat[K]: output matrix k \times m
 *              - const uint8_t rho[]: byte array containing seed rho
 **************************************************/
void polymatkm_expand(polyvecm mat[K], const uint8_t rho[SEEDBYTES]) {
    unsigned int i, j;

    for (i = 0; i < K; ++i)
        for (j = 0; j < M; ++j)
            poly_uniform(&mat[i].vec[j], rho, (i << 8) + j);
}

// doubles k * m sub-matrix of k * l mat
void polymatkl_double(polyvecl mat[K]) {
    unsigned int i, j, k;
    for (i = 0; i < K; ++i) {
        for (j = 1; j < L; ++j) {
            for (k = 0; k < N; ++k) {
                mat[i].vec[j].coeffs[k] *= 2;
            }
        }
    }
}

void polymatkl_pointwise_montgomery(polyveck *t, const polyvecl mat[K],
                                    const polyvecl *v) {
    unsigned int i;

    for (i = 0; i < K; ++i) {
        polyvecl_pointwise_acc_montgomery(&t->vec[i], &mat[i], v);
    }
}

void polymatkm_pointwise_montgomery(polyveck *t, const polyvecm mat[K],
                                    const polyvecm *v) {
    unsigned int i;

    for (i = 0; i < K; ++i) {
        polyvecm_pointwise_acc_montgomery(&t->vec[i], &mat[i], v);
    }
}


/*
 * Streamed version of polymatkl_pointwise_montgomery:
 * computes t = A1 * v where A1 = (col0 || 2*A0) and A0 entries are generated
 * on the fly from rho via poly_uniform(...).
 *
 * - col0 must be in NTT domain already.
 * - v must be in NTT domain.
 * - output t is in NTT domain.
 */
void polymatkl_pointwise_montgomery_stream(polyveck *t, const uint8_t rho[SEEDBYTES],
                                          const polyveck *col0, const polyvecl *v) {
    unsigned int i, j, k;

#if !FROZEN_A
    poly a;
#else
    poly_frozen af;
#endif

    for (i = 0; i < K; ++i) {
        /* column 0 */
        poly_pointwise_montgomery(&t->vec[i], &col0->vec[i], &v->vec[0]);

        /* columns 1..M (which correspond to 2*A0) */
        for (j = 1; j < L; ++j) {
            const uint16_t nonce = (uint16_t)((i << 8) + (j - 1));
#if !FROZEN_A
            poly_uniform(&a, rho, nonce);
            for (k = 0; k < N; ++k)
                a.coeffs[k] *= 2;
            poly_pointwise_montgomery_addto(&t->vec[i], &a, &v->vec[j]);
#else
            poly_uniform_frozen(&af, rho, nonce);
            /* multiply by 2 modulo q to keep uint16_t range */
            for (k = 0; k < N; ++k) {
                uint32_t x = (uint32_t)af.coeffs[k] * 2u;
                if (x >= Q) x -= Q;
                af.coeffs[k] = (uint16_t)x;
            }
            poly_pointwise_montgomery_mixed_addto(&t->vec[i], &af, &v->vec[j]);
#endif
        }
    }
}
