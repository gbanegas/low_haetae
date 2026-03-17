#include "keygen.h"
#include "reduce.h"   /* MONT, MONTSQ, montgomery_reduce, caddq */

/*************************************************
 * Name:        crypto_sign_keypair
 *
 * Description: Generates public and private key.
 *
 * Arguments:   - uint8_t *pk: pointer to output public key (allocated
 *                             array of CRYPTO_PUBLICKEYBYTES bytes)
 *              - uint8_t *sk: pointer to output private key (allocated
 *                             array of CRYPTO_SECRETKEYBYTES bytes)
 *
 * Returns 0 (success)
 **************************************************/

int crypto_sign_keypair(uint8_t *pk, uint8_t *sk)
{
    uint8_t seedbuf[2 * SEEDBYTES + CRHBYTES] = {0};
    uint16_t counter = 0;
    const uint8_t *rhoprime, *sigma, *key;
    polyvecm A[K], s1, s1hat;
    polyveck b, s2;
#if D > 0
    polyveck a, b0;
#else
    polyveck s2hat;
#endif

    // Get entropy \rho
    randombytes(seedbuf, SEEDBYTES);

    generate_seed_from_one_source(seedbuf, 2 * SEEDBYTES + CRHBYTES, seedbuf,
	SEEDBYTES);

    rhoprime = seedbuf;
    sigma = rhoprime + SEEDBYTES;
    key = sigma + CRHBYTES;

    // Expand Matrix A0 and vector a
    polymatkm_expand(A, rhoprime);


#ifdef ENABLE_KEYPAIR_MATRIX_BUFFER
    // Expand Matrix A0 and vector a
    for (size_t row = 0; row < K; ++row){
        for (size_t column = 0; column < M; ++column){
            poly_uniform_frozen(&A_gen[row].vec[column], rhoprime, nonce_A_gen(row, column));
        }
    }
#endif
#if D > 0
    /**********************************************
     * If there is rounding (D > 0), we need another polyveck a.
     * Then, b = a + A0 * s1 + s2 and the lower D bits are
     * rounded from b. The lower D bits are subsequently
     * subtracted from s2.
     **********************************************/
    polyveck_expand(&a, rhoprime);

reject:
    // Sample secret vectors s1 and s2
    polyvecmk_uniform_eta(&s1, &s2, sigma, counter);
    counter += M + K;

    // b = a + A0 * s1 + s2 mod q
    s1hat = s1;
    polyvecm_ntt(&s1hat);
    polymatkm_pointwise_montgomery(&b, A, &s1hat);
    polyveck_invntt_tomont(&b);
    polyveck_add(&b, &b, &s2);
    polyveck_add(&b, &b, &a);
    polyveck_freeze(&b);

    // round off D bits
    polyveck_decompose_vk(&b0, &b);
    polyveck_sub(&s2, &s2, &b0);

    int64_t squared_singular_value = polyvecmk_sqsing_value(&s1, &s2);
    if (squared_singular_value > GAMMA * GAMMA * N)
    {
        goto reject;
    }
#else
    /**********************************************
     * If there is no rounding (D == 0), we store
     * -2b directly in NTT domain into the public key.
     **********************************************/
reject:
    // Sample secret vectors s1 and s2
    polyvecmk_uniform_eta(&s1, &s2, sigma, counter);
    counter += M + K;
    int64_t squared_singular_value = polyvecmk_sqsing_value(&s1, &s2);
    if (squared_singular_value > GAMMA * GAMMA * N)
    {
        goto reject;
    }

    // b = A0 * s1 + s2 mod q
    s1hat = s1;
    s2hat = s2;
    polyvecm_ntt(&s1hat);
    polyveck_ntt(&s2hat);
    polymatkm_pointwise_montgomery(&b, A, &s1hat);
    /*
     * In D=0, public keys store the first column of A1 in NTT/Montgomery form.
     * Here, A0 coefficients are sampled in standard representation, while
     * s1hat/s2hat are in NTT+Montgomery. The pointwise multiplication therefore
     * yields standard-domain values; we must convert to Montgomery before
     * adding s2hat (and before packing into the public key).
     */
    polyveck_frommont(&b);
    polyveck_add(&b, &b, &s2hat);
    polyveck_double_negate(&b);
    polyveck_caddq(&b); // directly compute and store NTT(-2b)
#endif

    pack_pk(pk, &b, rhoprime);
    pack_sk(sk, pk, &s1, &s2, key);

    return 0;
}


static inline uint16_t nonce_a(const size_t row) {
	return (K << 8) + M + row;
}
static inline uint16_t nonce_A_gen(const size_t row, const size_t column) {
	return (row << 8) + column;
}
static inline uint16_t nonce_s_gen(const size_t nonce_offset,
		const size_t column) {
	return nonce_offset + column;
}
static inline uint16_t nonce_e_gen(const size_t nonce_offset, const size_t row) {
	return nonce_offset + M + row;
}

#if FROZEN_A == 1
/* Reusable scratch for streamed A_{row,col} expansion.
 * Keeping it at file scope avoids putting 512 bytes on the stack.
 * (This code path is single-threaded on Cortex-M4/QEMU.)
 */
static poly_frozen g_A_frozen_tmp;
#endif
/* FROZEN_A >= 2: fused sample+multiply eliminates g_A_frozen_tmp entirely. */

static void pointwise_montgomery_by_A_gen_elem_frozen_set(poly *dest,
		const uMatrixPointerM_frozen agenptr, const size_t row,
		const size_t column, const poly *src_hat) {
#ifdef ENABLE_KEYPAIR_MATRIX_BUFFER
    poly_pointwise_montgomery_mixed(dest, &agenptr.vec[row].vec[column], src_hat);
#else
#if FROZEN_A >= 2
    poly_uniform_frozen_pointwise_montgomery_set(dest, agenptr.seed, nonce_A_gen(row, column), src_hat);
#elif FROZEN_A == 1
    poly_uniform_frozen(&g_A_frozen_tmp, agenptr.seed, nonce_A_gen(row, column));
    poly_pointwise_montgomery_mixed(dest, &g_A_frozen_tmp, src_hat);
#else
    (void)agenptr; (void)row; (void)column; (void)src_hat;
#endif
#endif /* ENABLE_KEYPAIR_MATRIX_BUFFER */
}

static void pointwise_montgomery_by_A_gen_elem_frozen_addto(poly *acc,
		const uMatrixPointerM_frozen agenptr, const size_t row,
		const size_t column, const poly *src_hat) {
#ifdef ENABLE_KEYPAIR_MATRIX_BUFFER
    poly_pointwise_montgomery_mixed_addto(acc, &agenptr.vec[row].vec[column], src_hat);
#else
#if FROZEN_A >= 2
    poly_uniform_frozen_pointwise_montgomery_addto(acc, agenptr.seed, nonce_A_gen(row, column), src_hat);
#elif FROZEN_A == 1
    poly_uniform_frozen(&g_A_frozen_tmp, agenptr.seed, nonce_A_gen(row, column));
    poly_pointwise_montgomery_mixed_addto(acc, &g_A_frozen_tmp, src_hat);
#else
    (void)agenptr; (void)row; (void)column; (void)src_hat;
#endif
#endif /* ENABLE_KEYPAIR_MATRIX_BUFFER */
}


/*************************************************
 * Name:        crypto_sign_keypair
 *
 * Description: Generates public and private key.
 *
 * Arguments:   - uint8_t *pk: pointer to output public key (allocated
 *                             array of CRYPTO_PUBLICKEYBYTES bytes)
 *              - uint8_t *sk: pointer to output private key (allocated
 *                             array of CRYPTO_SECRETKEYBYTES bytes)
 *
 * Returns 0 (success)
 **************************************************/

int crypto_sign_keypair_low_stack(uint8_t *pk, uint8_t *sk)
{
    uint8_t seedbuf[2 * SEEDBYTES + CRHBYTES] = {0};
    uint16_t counter = 0;
    const uint8_t *rhoprime, *sigma, *key;
    polyvecm s1;
    polyveck b, s2;
    uMatrixPointerM_frozen A_gen_ptr = { .seed = seedbuf };
    /* scratch reused across NTT, adding 'a', and extracting b0 */
    poly tmp;
#if D > 0
    /* no additional large temporaries: 'a' and b0 are streamed per-row */
#else
    polyveck s2hat;
    /* D==0: stream A0*s1hat without materializing A on the stack. */
#if !FROZEN_A
    /* When not using frozen-A, we need a 1KB scratch poly for streamed A_{row,col}. */
    poly a;
#endif
#endif

    // Get entropy \rho
    randombytes(seedbuf, SEEDBYTES);

    generate_seed_from_one_source(seedbuf, 2 * SEEDBYTES + CRHBYTES, seedbuf,
	SEEDBYTES);

    rhoprime = seedbuf;
    sigma = rhoprime + SEEDBYTES;
    key = sigma + CRHBYTES;

    // Expand Matrix A0: streamed on-the-fly for both D>0 and D==0.

#if D > 0
    /**********************************************
     * If there is rounding (D > 0), we need another polyveck a.
     * Then, b = a + A0 * s1 + s2 and the lower D bits are
     * rounded from b. The lower D bits are subsequently
     * subtracted from s2.
     **********************************************/
    /* 'a' is sampled per-row later to avoid keeping a full polyveck on stack. */
    //sp_scope_t sc;

reject:
    // Sample secret vectors s1 and s2

   // polyvecmk_uniform_eta(&s1, &s2, sigma, counter);

    // Sample secret vectors s_gen and e_gen
    //sc = sp_begin("keypair:poly_uniform_eta");
	for (size_t column = 0; column < M; column++) {
		poly_uniform_eta(&s1.vec[column], sigma, nonce_s_gen(counter, column));
	}
    //sp_end(&sc);

	for (size_t row = 0; row < K; row++) {
		poly_uniform_eta(&s2.vec[row], sigma, nonce_e_gen(counter, row));
	}
  
    counter += M + K;

    /* b = a + A0*s1 + s2 (mod q), streaming A and avoiding s1hat/temp polys. */
    for (size_t column = 0; column < M; column++) {
        tmp = s1.vec[column];
       // sc = sp_begin("keypair:poly_ntt");
        poly_ntt(&tmp);
       // sp_end(&sc);

       // sc = sp_begin("keypair:pointwise_montgomery_by_A_gen_elem_frozen");
        for (size_t row = 0; row < K; row++) {
            if (column == 0) {
                pointwise_montgomery_by_A_gen_elem_frozen_set(&b.vec[row], A_gen_ptr,
                        row, column, &tmp);
            } else {
                pointwise_montgomery_by_A_gen_elem_frozen_addto(&b.vec[row], A_gen_ptr,
                        row, column, &tmp);
            }
        }
       // sp_end(&sc);
    }
    for (size_t row = 0; row < K; row++) {
		poly *b_elem = &b.vec[row];
		poly *e_gen_elem = &s2.vec[row];

		poly_invntt_tomont(b_elem);
		poly_add(b_elem, b_elem, e_gen_elem);

		/* add a-row (streamed) */
		poly_uniform(&tmp, rhoprime, nonce_a(row));
		poly_add(b_elem, b_elem, &tmp);
		poly_freeze(b_elem);

		/* round off D bits: b0 <- low bits of b, then s2 -= b0 */
		poly_decompose_vk(&tmp, b_elem);
		poly_sub(e_gen_elem, e_gen_elem, &tmp);
	}

	int64_t squared_singular_value = polyvecmk_sqsing_value(&s1, &s2);
	if (squared_singular_value > GAMMA * GAMMA * N) {
		goto reject;
	}
   // polymatkm_pointwise_montgomery(&b, A, &s1hat);
    
   /* polyveck_invntt_tomont(&b);
   
    polyveck_add(&b, &b, &s2);
    
    polyveck_add(&b, &b, &a);
    
    polyveck_freeze(&b);
  

    // round off D bits
    polyveck_decompose_vk(&b0, &b);
    polyveck_sub(&s2, &s2, &b0);

    int64_t squared_singular_value = polyvecmk_sqsing_value(&s1, &s2);
    if (squared_singular_value > GAMMA * GAMMA * N)
    {
        goto reject;
    }
        */
#else
    /**********************************************
     * If there is no rounding (D == 0), we store
     * -2b directly in NTT domain into the public key.
     **********************************************/
reject:
    // Sample secret vectors s1 and s2
    polyvecmk_uniform_eta(&s1, &s2, sigma, counter);
    counter += M + K;
    int64_t squared_singular_value = polyvecmk_sqsing_value(&s1, &s2);
    if (squared_singular_value > GAMMA * GAMMA * N)
    {
        goto reject;
    }

    // b = A0 * s1 + s2 mod q, then store NTT(-2b) in pk.
    //
    // We keep the same representation conventions as the reference:
    // - s1 is sampled in coefficient form, then we NTT each column on demand.
    // - A0 entries are generated on-the-fly (either poly_uniform or poly_uniform_frozen)
    //   and interpreted directly as NTT-domain polynomials (NTT is a bijection).
    // - s2 is NTT'd into s2hat.
    s2hat = s2;
    polyveck_ntt(&s2hat);

    for (size_t column = 0; column < M; ++column) {
        tmp = s1.vec[column];
        poly_ntt(&tmp);

        for (size_t row = 0; row < K; ++row) {
            const uint16_t nonce = nonce_A_gen(row, column);
#if FROZEN_A >= 2
            if (column == 0) {
                poly_uniform_frozen_pointwise_montgomery_set(&b.vec[row], A_gen_ptr.seed, nonce, &tmp);
            } else {
                poly_uniform_frozen_pointwise_montgomery_addto(&b.vec[row], A_gen_ptr.seed, nonce, &tmp);
            }
#elif FROZEN_A == 1
            poly_uniform_frozen(&g_A_frozen_tmp, A_gen_ptr.seed, nonce);
            if (column == 0) {
                poly_pointwise_montgomery_mixed(&b.vec[row], &g_A_frozen_tmp, &tmp);
            } else {
                poly_pointwise_montgomery_mixed_addto(&b.vec[row], &g_A_frozen_tmp, &tmp);
            }
#else
            poly_uniform(&a, rhoprime, nonce);
            if (column == 0) {
                poly_pointwise_montgomery(&b.vec[row], &a, &tmp);
            } else {
                poly_pointwise_montgomery_addto(&b.vec[row], &a, &tmp);
            }
#endif
        }
    }

    /* Convert to Montgomery domain before adding s2hat / packing. */
    polyveck_frommont(&b);
    polyveck_add(&b, &b, &s2hat);
    polyveck_double_negate(&b);
    polyveck_caddq(&b); // directly compute and store NTT(-2b)
#endif

    pack_pk(pk, &b, rhoprime);
    pack_sk(sk, pk, &s1, &s2, key);

    return 0;
}

/*
 * Streamed / low-stack key generation variants.
 *
 * These functions generate the same (pk,sk) format as crypto_sign_keypair_low_stack,
 * but reduce peak stack by:
 *   - sampling and packing s1/s2 polynomials on the fly (no polyvecm/polyveck secrets on stack),
 *   - streaming the random 'a' term directly into the accumulator (no temporary poly),
 *   - decomposing and packing b directly into pk (no extra b1/b0 temporaries).
 */

int crypto_sign_keypair_streamed_v1(uint8_t *pk, uint8_t *sk)
{
    uint8_t seedbuf[2 * SEEDBYTES + CRHBYTES] = {0};
    uint16_t counter = 0;
    const uint8_t *rhoprime, *sigma, *key;
    uMatrixPointerM_frozen A_gen_ptr = {.seed = seedbuf};

    polyveck b;
    poly s;
    int32_t sum[N];

#if D > 0
    /* no additional large temporaries: 'a' and b0 are streamed per-row */
#else
    /* D==0 path keeps the reference implementation for now. */
    polyveck s2hat;
    polyvecm A[K], s1hat;
    polyvecm s1;
    polyveck s2;
#endif

    /* Get entropy rho, then derive (rhoprime || sigma || key). */
    randombytes(seedbuf, SEEDBYTES);
    generate_seed_from_one_source(seedbuf, 2 * SEEDBYTES + CRHBYTES, seedbuf, SEEDBYTES);

    rhoprime = seedbuf;
    sigma = rhoprime + SEEDBYTES;
    key = sigma + CRHBYTES;

#if D > 0
reject:
    /* b = a + A0*s1 + s2 (mod q), streaming A and avoiding s1hat/temp polys. */
    memset(sum, 0, sizeof(int32_t) * N);

    for (size_t column = 0; column < M; column++)
    {
        poly_uniform_eta(&s, sigma, nonce_s_gen(counter, column));
        polyvecmk_sqsing_value_s(sum, &s);
        pack_sk_stream_s1(sk, s, column);
        poly_ntt(&s);

        for (size_t row = 0; row < K; row++)
        {
            if (column == 0)
            {
                pointwise_montgomery_by_A_gen_elem_frozen_set(&b.vec[row], A_gen_ptr, row, column, &s);
            }
            else
            {
                pointwise_montgomery_by_A_gen_elem_frozen_addto(&b.vec[row], A_gen_ptr, row, column, &s);
            }
        }
    }

    for (size_t row = 0; row < K; row++)
    {
        /* add s2-row (streamed) */
        poly_invntt_tomont(&b.vec[row]);
        poly_uniform_eta(&s, sigma, nonce_e_gen(counter, row));
        poly_add(&b.vec[row], &b.vec[row], &s);

        /* add a-row (streamed) */
        poly_uniform_add(&b.vec[row], rhoprime, nonce_a(row));
        poly_freeze(&b.vec[row]);

        /* round off D bits: b0 <- low bits of b, then s2 -= b0 */
        poly_decompose_stream_pack(pk, &b.vec[row], row);
        poly_sub(&s, &s, &b.vec[row]);

        polyvecmk_sqsing_value_s(sum, &s);
        pack_sk_stream_s2(sk, s, row);
    }

    int64_t squared_singular_value = polyvecmk_sqsing_value_min(sum, &s);

    counter += M + K;
    if (squared_singular_value > (int64_t)(GAMMA * GAMMA * N))
    {
        goto reject;
    }

    pack_pk_stream_remain(pk, rhoprime);
    pack_sk_stream_remain(sk, pk, key);
#else
    /**********************************************
     * If there is no rounding (D == 0), we store
     * -2b directly in NTT domain into the public key.
     **********************************************/
reject:
    /* Sample secret vectors s1 and s2 */
    polyvecmk_uniform_eta(&s1, &s2, sigma, counter);
    counter += M + K;
    int64_t squared_singular_value = polyvecmk_sqsing_value(&s1, &s2);
    if (squared_singular_value > GAMMA * GAMMA * N)
    {
        goto reject;
    }

    /* b = A0 * s1 + s2 mod q */
    s1hat = s1;
    s2hat = s2;
    polyvecm_ntt(&s1hat);
    polyveck_ntt(&s2hat);
    polymatkm_pointwise_montgomery(&b, A, &s1hat);
    polyveck_frommont(&b);
    polyveck_add(&b, &b, &s2hat);
    polyveck_double_negate(&b);
    polyveck_caddq(&b); /* directly compute and store NTT(-2b) */

    pack_pk(pk, &b, rhoprime);
    pack_sk(sk, pk, &s1, &s2, key);
#endif

    return 0;
}

int crypto_sign_keypair_streamed_v2(uint8_t *pk, uint8_t *sk)
{
    uint8_t seedbuf[2 * SEEDBYTES + CRHBYTES] = {0};
    uint16_t counter = 0;
    const uint8_t *rhoprime, *sigma, *key;
    uMatrixPointerM_frozen A_gen_ptr = {.seed = seedbuf};

    /*
     * NOTE:
     *  - For D>0 (rounded/truncated VK variants), we stream one row at a time, so we only
     *    need a single polynomial accumulator 'b' and a scratch 's', plus the singular-value
     *    accumulators.
     *  - For D==0 (no rounding), the reference path below computes the full vector b \in R^K,
     *    so 'b' must be a polyveck.
     */
    poly b, s;
    int32_t sum[N];

    /* Get entropy rho, then derive (rhoprime || sigma || key). */
    randombytes(seedbuf, SEEDBYTES);
    generate_seed_from_one_source(seedbuf, 2 * SEEDBYTES + CRHBYTES, seedbuf, SEEDBYTES);

    rhoprime = seedbuf;
    sigma = rhoprime + SEEDBYTES;
    key = sigma + CRHBYTES;

#if D > 0
reject:
    /* b = a + A0*s1 + s2 (mod q), streaming A and avoiding s1hat/temp polys. */
    memset(sum, 0, sizeof(int32_t) * N);

    for (size_t row = 0; row < K; row++)
    {
        for (size_t column = 0; column < M; column++)
        {
            /* re-sample s1[column] deterministically; pack s1 only once (row==0) */
            poly_uniform_eta(&s, sigma, nonce_s_gen(counter, column));
            if (row == 0)
            {
                polyvecmk_sqsing_value_s(sum, &s);
                pack_sk_stream_s1(sk, s, column);
            }
            poly_ntt(&s);

            if (column == 0)
            {
                pointwise_montgomery_by_A_gen_elem_frozen_set(&b, A_gen_ptr, row, column, &s);
            }
            else
            {
                pointwise_montgomery_by_A_gen_elem_frozen_addto(&b, A_gen_ptr, row, column, &s);
            }
        }

        /* add s2-row (streamed) */
        poly_invntt_tomont(&b);
        poly_uniform_eta(&s, sigma, nonce_e_gen(counter, row));
        poly_add(&b, &b, &s);

        /* add a-row (streamed) */
        poly_uniform_add(&b, rhoprime, nonce_a(row));
        poly_freeze(&b);

        /* round off D bits: b0 <- low bits of b, then s2 -= b0 */
        poly_decompose_stream_pack(pk, &b, row);
        poly_sub(&s, &s, &b);

        polyvecmk_sqsing_value_s(sum, &s);
        pack_sk_stream_s2(sk, s, row);
    }

    int64_t squared_singular_value = 0;
    poly scratch;
    squared_singular_value = polyvecmk_sqsing_value_min(sum, &scratch);

    counter += M + K;
    if (squared_singular_value > (int64_t)(GAMMA * GAMMA * N))
    {
        goto reject;
    }

    pack_pk_stream_remain(pk, rhoprime);
    pack_sk_stream_remain(sk, pk, key);
#else
	/**********************************************
	 * If there is no rounding (D == 0), we store
	 * -2b directly in NTT domain into the public key.
	 *
	 * IMPORTANT:
	 *   The previous block here relied on an uninitialized A[K] and therefore
	 *   produced invalid public keys (signatures never verify). We instead
	 *   compute b row-by-row while expanding A0 on the fly.
	 **********************************************/
reject:
	memset(sum, 0, sizeof(sum));
	uint16_t nonce_offset = counter;

	/* Sample and pack s1 (length M) */
	for (size_t column = 0; column < M; column++)
	{
		poly_uniform_eta(&s, sigma, nonce_s_gen(nonce_offset, column));
		polyvecmk_sqsing_value_s(sum, &s);
		pack_sk_stream_s1(sk, s, column);
	}

	/* Sample and pack s2 (length K) */
	for (size_t row = 0; row < K; row++)
	{
		poly_uniform_eta(&s, sigma, nonce_e_gen(nonce_offset, row));
		polyvecmk_sqsing_value_s(sum, &s);
		pack_sk_stream_s2(sk, s, row);
	}
    int64_t squared_singular_value = 0;
    poly scratch;
    squared_singular_value = polyvecmk_sqsing_value_min(sum, &scratch);

	/* Consume nonces for next attempt */
	counter += M + K;
	if (squared_singular_value > (int64_t)(GAMMA * GAMMA * N))
	{
		goto reject;
	}

	/* Compute and pack pk = (rhoprime, NTT(-2b)) row-by-row. */
	pack_pk_stream_remain(pk, rhoprime);
	for (size_t row = 0; row < K; row++)
	{
		/* b_hat = sum_j A[row,j] * NTT(s1[j])  (NTT domain, standard rep) */
		for (size_t column = 0; column < M; column++)
		{
			poly_uniform_eta(&s, sigma, nonce_s_gen(nonce_offset, column));
			poly_ntt(&s);

			if (column == 0)
			{
				pointwise_montgomery_by_A_gen_elem_frozen_set(&b, A_gen_ptr, row, column, &s);
			}
			else
			{
				pointwise_montgomery_by_A_gen_elem_frozen_addto(&b, A_gen_ptr, row, column, &s);
			}
		}

		/* Convert b_hat to Montgomery domain before adding s2hat. */
		for (size_t i = 0; i < N; i++)
		{
			b.coeffs[i] = montgomery_reduce((int64_t)b.coeffs[i] * MONTSQ);
		}

		/* Add s2hat(row) (NTT+Montgomery). */
		poly_uniform_eta(&s, sigma, nonce_e_gen(nonce_offset, row));
		poly_ntt(&s);
		poly_add(&b, &b, &s);

		/* Multiply by -2 in Montgomery domain and map into [0, Q-1]. */
		for (size_t i = 0; i < N; i++)
		{
			b.coeffs[i] = montgomery_reduce((int64_t)b.coeffs[i] * MONT * -2);
			b.coeffs[i] = caddq(b.coeffs[i]);
		}

		pack_pk_stream(pk, &b, row);
	}

	pack_sk_stream_remain(sk, pk, key);
#endif

    return 0;
}
