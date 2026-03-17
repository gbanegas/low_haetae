/*
 * poly_sparsemul.h
 *
 * Prototypes for low-stack sparse-challenge helpers.
 */
#ifndef POLY_SPARSEMUL_H
#define POLY_SPARSEMUL_H

#include <stdint.h>
#include "poly.h"
#include "polyfix.h"

#ifdef __cplusplus
extern "C" {
#endif

unsigned poly_c_to_sparse(uint16_t *idx, int8_t *sgn, unsigned max_w, const poly *c);

void poly_mul_sparse_negacyclic(poly *out, const poly *s,
                                const uint16_t *idx, const int8_t *sgn,
                                unsigned w);

uint64_t polyfix_addmul_sparse_inplace_acc2(polyfix *y, const poly *s,
                                           const uint16_t *idx, const int8_t *sgn,
                                           unsigned w,
                                           int32_t scale,
                                           uint64_t acc2, uint64_t bound0);

#ifdef __cplusplus
}
#endif

#endif /* POLY_SPARSEMUL_H */
