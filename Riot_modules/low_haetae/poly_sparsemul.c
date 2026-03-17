/*
 * poly_sparsemul.c
 *
 * Low-stack helpers for HAETAE signing:
 *   - Extract a sparse representation of the challenge polynomial c (dense poly -> (idx,sgn) list).
 *   - Compute (c * s) in R_q = Z_q[x]/(x^N + 1) using only negacyclic shifts (no NTT).
 *   - (Optional) Update a fixed-point polynomial y in-place with z = y + scale*(c*s),
 *     while simultaneously accumulating the "2z - y" squared-norm accumulator (acc2).
 *
 * Rationale:
 *   If c has small weight tau and coefficients in {+1,0,-1}, then
 *     c(x)s(x) = sum_{k=1..tau} sgn_k * x^{idx_k} * s(x)   (mod x^N+1)
 *   and multiplication by x^{idx_k} is a negacyclic rotation with a sign flip on wrap.
 *   Coefficients remain small (no wrap modulo q) when s has small coefficients.
 *
 * This file is pure C and intended for Cortex-M4 / low-stack builds.
 */

#include <stdint.h>
#include <stddef.h>

#include "params.h"
#include "poly.h"
#include "polyfix.h"

/* If your build uses namespacing via config.h, you can uncomment this:
 *
 * #include "config.h"
 * #define poly_c_to_sparse                    HAETAE_NAMESPACE(poly_c_to_sparse)
 * #define poly_mul_sparse_negacyclic          HAETAE_NAMESPACE(poly_mul_sparse_negacyclic)
 * #define polyfix_addmul_sparse_inplace_acc2  HAETAE_NAMESPACE(polyfix_addmul_sparse_inplace_acc2)
 */

/*
 * Normalize a challenge coefficient to {-1,0,+1}.
 * Many implementations store -1 as (Q-1) mod Q.
 */
static inline int8_t c_coeff_to_pm1(int32_t v) {
    if (v == 0) return 0;
    if (v == 1) return +1;
    if ((uint32_t)v == (uint32_t)(Q - 1)) return -1;
    /* If your challenge uses another encoding, adapt here. */
    return 0;
}

/*
 * Extract sparse list (idx[], sgn[]) from dense challenge polynomial c.
 *
 * - idx[] holds positions in [0, N-1]
 * - sgn[] holds +/-1 (int8_t)
 * - max_w is the maximum number of nonzero entries to record (typically TAU)
 *
 * Returns:
 *   w = number of entries recorded (<= max_w).
 *
 * NOTE:
 *   If c is guaranteed to have exactly TAU nonzeros, w should equal TAU.
 */
unsigned poly_c_to_sparse(uint16_t *idx, int8_t *sgn, unsigned max_w, const poly *c) {
    unsigned w = 0;
    for (unsigned i = 0; i < N; i++) {
        int8_t s = c_coeff_to_pm1(c->coeffs[i]);
        if (!s) continue;

        idx[w] = (uint16_t)i;
        sgn[w] = s;
        if (++w == max_w) break;
    }
    return w;
}

/*
 * out = c*s in R_q = Z_q[x]/(x^N + 1), using sparse (idx,sgn) list for c.
 *
 * This produces SIGNED coefficients (small range) and performs no reduction mod q.
 * This is appropriate when:
 *   - s is small (eta-distribution),
 *   - c is sparse +/-1,
 *   - and you will use the product for fixed-point additions / norm checks
 *     where a centered small representation is desired.
 *
 * If you need canonical [0,Q) representatives, add a reduction step afterward.
 */
void poly_mul_sparse_negacyclic(poly *out, const poly *s,
                                const uint16_t *idx, const int8_t *sgn,
                                unsigned w) {
    for (unsigned i = 0; i < N; i++) out->coeffs[i] = 0;

    for (unsigned k = 0; k < w; k++) {
        const unsigned sh = (unsigned)idx[k];
        const int32_t cs = (int32_t)sgn[k]; /* +/-1 */

        /* Term: cs * x^sh * s(x)  (mod x^N+1) */
        for (unsigned j = 0; j < N; j++) {
            unsigned dst = j + sh;
            if (dst < N) {
                out->coeffs[dst] += cs * s->coeffs[j];
            } else {
                dst -= N;
                /* wrap => multiply by x^N = -1 */
                out->coeffs[dst] -= cs * s->coeffs[j];
            }
        }
    }
}

/*
 * In-place low-stack update:
 *   y <- y + scale * (c*s)
 * and optionally accumulate:
 *   u = y_old + 2*scale*(c*s)  (i.e., 2z - y_old)
 *   acc2 <- min(bound0, acc2 + sum_i u_i^2)
 *
 * Arguments:
 *   - y:      fixed-point poly (updated in place)
 *   - s:      small secret poly (eta)
 *   - idx,sgn,w: sparse challenge representation
 *   - scale:  typically (sgn_from_b * LN)
 *   - acc2:   current accumulator (caller maintains it)
 *   - bound0: saturation threshold (e.g., B0SQ*LN^2)
 *
 * Returns:
 *   updated acc2.
 *
 * Notes:
 *   - This avoids allocating a temporary poly "tmp" on the stack.
 *   - Complexity: O(w*N). For small w (tau) this is acceptable.
 */
uint64_t polyfix_addmul_sparse_inplace_acc2(polyfix *y, const poly *s,
                                           const uint16_t *idx, const int8_t *sgn,
                                           unsigned w,
                                           int32_t scale,
                                           uint64_t acc2, uint64_t bound0) {
    for (unsigned i = 0; i < N; i++) {
        int32_t sum = 0;

        /* Compute coefficient i of (c*s) in negacyclic ring. */
        for (unsigned k = 0; k < w; k++) {
            const unsigned sh = (unsigned)idx[k];
            const int32_t ck = (int32_t)sgn[k]; /* +/-1 */

            /* coefficient i of x^sh*s:
             * if i >= sh:  + s[i - sh]
             * else:       - s[i - sh + N]
             */
            if (i >= sh) {
                sum += ck * s->coeffs[i - sh];
            } else {
                sum -= ck * s->coeffs[i - sh + N];
            }
        }

        /* t = scale * (c*s)[i] */
        const int32_t y0 = y->coeffs[i];
        const int32_t t  = (int32_t)((int64_t)scale * (int64_t)sum);
        const int32_t z  = y0 + t;

        /* u = 2z - y = y + 2t */
        const int32_t u  = y0 + (t << 1);

        if (acc2 < bound0) {
            const uint64_t term2 = (uint64_t)((int64_t)u * (int64_t)u);
            if (term2 >= bound0 - acc2) acc2 = bound0;
            else acc2 += term2;
        }

        y->coeffs[i] = z;
    }

    return acc2;
}
