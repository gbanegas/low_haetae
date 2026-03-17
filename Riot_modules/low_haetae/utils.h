#ifndef HAETAE_UTILS_H_
#define HAETAE_UTILS_H_

#include "config.h"
#include "params.h"
#include "poly.h"
#include "polyfix.h"
#include "polyvec.h"
#include "polymat.h"
#include "fips202.h"
#include "packing.h"
#include "symmetric.h"

typedef union MatrixPointerL_frozen{
    const uint8_t *seed;
    polyvecl_frozen *vec;
} uMatrixPointerL_frozen;

typedef union MatrixPointerM_frozen{
    const uint8_t *seed;
    polyvecm_frozen *vec;
} uMatrixPointerM_frozen;


void generate_seed_from_one_source(uint8_t *seed, size_t seed_len,
		const uint8_t *src, size_t src_len) ;

#endif