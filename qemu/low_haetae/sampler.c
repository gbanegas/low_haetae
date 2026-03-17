#include "sampler.h"
#include "fixpoint.h"
#include "symmetric.h"
#include <stdint.h>
#include <string.h>

/*************************************************
 * Name:        rej_uniform
 *
 * Description: Sample uniformly random coefficients in [0, Q-1] by
 *              performing rejection sampling on array of random bytes.
 *
 * Arguments:   - int32_t *a: pointer to output array (allocated)
 *              - unsigned int len: number of coefficients to be sampled
 *              - const uint8_t *buf: array of random bytes
 *              - unsigned int buflen: length of array of random bytes
 *
 * Returns number of sampled coefficients. Can be smaller than len if not enough
 * random bytes were given.
 **************************************************/
unsigned int rej_uniform(int32_t *a, unsigned int len, const uint8_t *buf,
                         unsigned int buflen) {
    unsigned int ctr, pos;
    uint32_t t;

    ctr = pos = 0;
    while (ctr < len && pos + 2 <= buflen) {
        t = buf[pos++];
        t |= (uint32_t)buf[pos++] << 8;

        if (t < Q)
            a[ctr++] = t;
    }
    return ctr;
}
unsigned int rej_uniform_add(int32_t *a, unsigned int len, const uint8_t *buf,
                         unsigned int buflen) {
    unsigned int ctr, pos;
    uint32_t t;

    ctr = pos = 0;
    while (ctr < len && pos + 2 <= buflen) {
        t = buf[pos++];
        t |= (uint32_t)buf[pos++] << 8;

        if (t < Q)
            a[ctr++] += t;
    }
    return ctr;
}


#if FROZEN_A >= 1
unsigned int rej_uniform_frozen(uint16_t *a, unsigned int len, const uint8_t *buf,
                                unsigned int buflen) {
    unsigned int ctr, pos;
    uint32_t t;

    ctr = pos = 0;
    while (ctr < len && pos + 2 <= buflen) {
        t = buf[pos++];
        t |= (uint32_t)buf[pos++] << 8;
        if (t < Q)
            a[ctr++] = (uint16_t)t;
    }
    return ctr;
}
#endif



/*************************************************
 * Name:        rej_eta
 *
 * Description: Sample uniformly random coefficients in [-ETA, ETA] by
 *              performing rejection sampling on array of random bytes.
 *
 * Arguments:   - int32_t *a: pointer to output array (allocated)
 *              - unsigned int len: number of coefficients to be sampled
 *              - const uint8_t *buf: array of random bytes
 *              - unsigned int buflen: length of array of random bytes
 *
 * Returns number of sampled coefficients. Can be smaller than len if not enough
 * random bytes were given.
 **************************************************/
static int32_t mod3(uint8_t t) {
    int32_t r;
    r = (t >> 4) + (t & 0xf);
    r = (r >> 2) + (r & 3);
    r = (r >> 2) + (r & 3);
    r = (r >> 2) + (r & 3);
    return r - (3 * (r >> 1));
}
static int32_t mod3_leq26(uint8_t t) {
    int32_t r;
    r = (t >> 4) + (t & 0xf);
    r = (r >> 2) + (r & 3);
    r = (r >> 2) + (r & 3);
    return r - (3 * (r >> 1));
}
static int32_t mod3_leq8(uint8_t t) {
    int32_t r;
    r = (t >> 2) + (t & 3);
    r = (r >> 2) + (r & 3);
    return r - (3 * (r >> 1));
}
unsigned int rej_eta(int32_t *a, unsigned int len, const uint8_t *buf,
                     unsigned int buflen) {
    unsigned int ctr, pos;

    ctr = pos = 0;
    while (ctr < len && pos < buflen) {
#if ETA == 1
        uint32_t t = buf[pos++];
        if (t < 243) {
            // reduce mod 3
            a[ctr++] = mod3(t);

            if (ctr >= len)
                break;

            t *= 171; // 171*3 = 1 mod 256
            t >>= 9;
            a[ctr++] = mod3(t);

            if (ctr >= len)
                break;

            t *= 171;
            t >>= 9;
            a[ctr++] = mod3_leq26(t);

            if (ctr >= len)
                break;

            t *= 171;
            t >>= 9;
            a[ctr++] = mod3_leq8(t);

            if (ctr >= len)
                break;

            t *= 171;
            t >>= 9;
            a[ctr++] = (int32_t)t - (int32_t)3 * (t >> 1);
        }
#elif ETA == 2
        uint32_t t0, t1;
        t0 = buf[pos] & 0x0F;
        t1 = buf[pos++] >> 4;
        if (t0 < 15) {
            t0 = t0 - (205 * t0 >> 10) * 5;
            a[ctr++] = 2 - t0;
        }
        if (t1 < 15 && ctr < len) {
            t1 = t1 - (205 * t1 >> 10) * 5;
            a[ctr++] = 2 - t1;
        }
#endif
    }
    return ctr;
}


unsigned int rej_eta_add(int32_t *a, unsigned int len, const uint8_t *buf,
                     unsigned int buflen) {
    unsigned int ctr, pos;

    ctr = pos = 0;
    while (ctr < len && pos < buflen) {
#if ETA == 1
        uint32_t t = buf[pos++];
        if (t < 243) 
        {
            // reduce mod 3
            a[ctr++] += mod3(t);

            if (ctr >= len)
                break;

            t *= 171; // 171*3 = 1 mod 256
            t >>= 9;
            a[ctr++] += mod3(t);

            if (ctr >= len)
                break;

            t *= 171;
            t >>= 9;
            a[ctr++] += mod3_leq26(t);

            if (ctr >= len)
                break;

            t *= 171;
            t >>= 9;
            a[ctr++] += mod3_leq8(t);

            if (ctr >= len)
                break;

            t *= 171;
            t >>= 9;
            a[ctr++] += (int32_t)t - (int32_t)3 * (t >> 1);
        }
#elif ETA == 2
        uint32_t t0, t1;
        t0 = buf[pos] & 0x0F;
        t1 = buf[pos++] >> 4;
        if (t0 < 15) {
            t0 = t0 - (205 * t0 >> 10) * 5;
            a[ctr++] += 2 - t0;
        }
        if (t1 < 15 && ctr < len) {
            t1 = t1 - (205 * t1 >> 10) * 5;
            a[ctr++] += 2 - t1;
        }
#endif
    }
    return ctr;
}

static uint64_t approx_exp(const uint64_t x) {
    int64_t result;
    result = -0x0000B6C6340925AELL;
    result = ((smulh48(result, x) + (1LL << 2)) >> 3) + 0x0000B4BD4DF85227LL;
    result = ((smulh48(result, x) + (1LL << 2)) >> 3) - 0x0000887F727491E2LL;
    result = ((smulh48(result, x) + (1LL << 1)) >> 2) + 0x0000AAAA643C7E8DLL;
    result = ((smulh48(result, x) + (1LL << 1)) >> 2) - 0x0000AAAAA98179E6LL;
    result = ((smulh48(result, x) + 1LL) >> 1) + 0x0000FFFFFFFB2E7ALL;
    result = ((smulh48(result, x) + 1LL) >> 1) - 0x0000FFFFFFFFF85FLL;
    result = ((smulh48(result, x))) + 0x0000FFFFFFFFFFFCLL;
    return result;
}

#define CDTLEN 64
static const uint32_t CDT[CDTLEN] = {
 3266,  6520,  9748, 12938, 16079, 19159, 22168, 25096,
27934, 30674, 33309, 35833, 38241, 40531, 42698, 44742,
46663, 48460, 50135, 51690, 53128, 54454, 55670, 56781,
57794, 58712, 59541, 60287, 60956, 61554, 62085, 62556,
62972, 63337, 63657, 63936, 64178, 64388, 64569, 64724,
64857, 64970, 65066, 65148, 65216, 65273, 65321, 65361,
65394, 65422, 65444, 65463, 65478, 65490, 65500, 65508,
65514, 65519, 65523, 65527, 65529, 65531, 65533, 65534
};


static uint64_t sample_gauss16(const uint64_t rand16) {
    unsigned int i;
    uint64_t r = 0;
    for (i = 0; i < CDTLEN; i++) {
        r += (((uint64_t)CDT[i] - rand16) >> 63) & 1;
    }
    return r;
}

#define GAUSS_RAND (72 + 16 + 48)
#define GAUSS_RAND_BYTES ((GAUSS_RAND + 7) / 8)
static int sample_gauss_sigma76(uint64_t *r, fp96_76 *sqr,
                                const uint8_t rand[GAUSS_RAND_BYTES]) {
    const uint64_t rand_gauss16 = rand[0] | (((uint64_t) rand[1]) << 8); 
    const uint64_t rand_rej = rand[2] | (((uint64_t) rand[3]) << 8) | (((uint64_t) rand[4]) << 16) | (((uint64_t) rand[5]) << 24)
     | (((uint64_t) rand[6]) << 32) | (((uint64_t) rand[7]) << 40);
    uint64_t x, exp_in;
    fp96_76 y;

    // sample x
    x = sample_gauss16(rand_gauss16);

    // y := append x to y
    // leave 16 bit for carries
    y.limb48[0] = rand[8] | ((uint64_t)rand[9] << 8) |
                  ((uint64_t)rand[10] << 16) | ((uint64_t)rand[11] << 24) |
                  ((uint64_t)rand[12] << 32) | ((uint64_t)rand[13] << 40);
    y.limb48[1] =
        rand[14] | ((uint64_t)rand[15] << 8) | ((uint64_t)rand[16] << 16) |
        (x << 24);

    // r := round y 
    *r = (y.limb48[0] >> 15) ^ (y.limb48[1] << 33);
    *r += 1; // rounding
    *r >>= 1;

    // sqr := y*y
    fixpoint_square(sqr, &y);

    // sqr[1] = y^2 >> (76+48)                // 34 bit
    // sqr[0] = (y^2 >> 76) & ((1UL<<48)-1)   // 48 bit
    // exp_in := sqr - ((x*x) << 68)
    exp_in = sqr->limb48[1] - ((x*x) << (68 - 48));
    exp_in <<= 20;
    exp_in |= sqr->limb48[0] >> 28;
    exp_in += 1; // rounding
    exp_in >>= 1;

    return ((((int64_t)(rand_rej ^
                        (rand_rej & 1)) // set lowest bit to zero in order to
                                        // use it for rejection if sample==0
              - (int64_t)approx_exp(exp_in)) >>
             63) // reject with prob 1-approx_exp(exp_in)
            & (((*r | -*r) >> 63) | rand_rej)) &
           1; // if the sample is zero, clear the return value with prob 1/2
}

int sample_gauss(uint64_t *r, fp96_76 *sqsum, const uint8_t *buf, const size_t buflen, const size_t len, const int dont_write_last)
{
    const uint8_t *pos = buf;
    fp96_76 sqr;
    size_t bytecnt = buflen, coefcnt = 0, cnt = 0;
    int accepted;
    uint64_t dummy;
    
    while (coefcnt < len) {
        if (bytecnt < GAUSS_RAND_BYTES) {
          renormalize(sqsum);
          return coefcnt;
        }

        if (dont_write_last && coefcnt == len-1)
        {
          accepted = sample_gauss_sigma76(&dummy, &sqr, pos);
        } else {
          accepted = sample_gauss_sigma76(&r[coefcnt], &sqr, pos);
        }
        cnt += 1;
        coefcnt += accepted;
        pos += GAUSS_RAND_BYTES;
        bytecnt -= GAUSS_RAND_BYTES;

        sqsum->limb48[0] += sqr.limb48[0] & -(int64_t)accepted;
        sqsum->limb48[1] += sqr.limb48[1] & -(int64_t)accepted;
    }

    renormalize(sqsum);
    return len;
}

/*
 * Byte-wise squeeze wrapper around stream256_squeezeblocks().
 * This keeps portability across SHAKE-based and AES-CTR based streams,
 * while allowing us to avoid large stack buffers.
 */
typedef struct {
    stream256_state st;
    uint8_t block[STREAM256_BLOCKBYTES];
    size_t idx;
} stream256_bytectx;

static void stream256_bytectx_init(stream256_bytectx *ctx,
                                  const uint8_t seed[CRHBYTES],
                                  uint16_t nonce) {
    stream256_init(&ctx->st, seed, nonce);
    ctx->idx = STREAM256_BLOCKBYTES; /* mark as empty */
}

static void stream256_squeeze_bytes(stream256_bytectx *ctx,
                                    uint8_t *out, size_t outlen) {
    while (outlen > 0) {
        if (ctx->idx >= STREAM256_BLOCKBYTES) {
            stream256_squeezeblocks(ctx->block, 1, &ctx->st);
            ctx->idx = 0;
        }
        size_t take = STREAM256_BLOCKBYTES - ctx->idx;
        if (take > outlen)
            take = outlen;
        memcpy(out, ctx->block + ctx->idx, take);
        ctx->idx += take;
        out += take;
        outlen -= take;
    }
}

/*
 * Low-stack streaming sampler for N Gaussian magnitudes.
 *
 * We keep a small buffer and only carry the < GAUSS_RAND_BYTES remainder
 * between refills, preserving the exact byte stream order.
 */
void sample_gauss_N(uint64_t *r, uint8_t *signs, fp96_76 *sqsum,
                    const uint8_t seed[CRHBYTES], const uint16_t nonce,
                    const size_t len) {
    stream256_bytectx ctx;
    stream256_bytectx_init(&ctx, seed, nonce);

    const size_t sign_bytes = len / 8;
    if (sign_bytes > 0)
        stream256_squeeze_bytes(&ctx, signs, sign_bytes);

    /* Small refill buffer. */
#define GAUSS_BUF_BYTES (STREAM256_BLOCKBYTES + GAUSS_RAND_BYTES)
    uint8_t buf[GAUSS_BUF_BYTES];
    size_t bytecnt = 0;
    size_t coefcnt = 0;
    const int dont_write_last = (len % N);

    while (coefcnt < len) {
        /* Refill buffer if we cannot form a complete GAUSS_RAND_BYTES chunk. */
        if (bytecnt < GAUSS_RAND_BYTES) {
            /* Keep the remainder (bytecnt < GAUSS_RAND_BYTES) at the front. */
            const size_t need = GAUSS_BUF_BYTES - bytecnt;
            stream256_squeeze_bytes(&ctx, buf + bytecnt, need);
            bytecnt += need;
        }

        /* Consume as many full GAUSS_RAND_BYTES chunks as available. */
        coefcnt += (size_t)sample_gauss(r + coefcnt, sqsum, buf, bytecnt,
                                        len - coefcnt, dont_write_last);

        /* Preserve only the remainder (< GAUSS_RAND_BYTES) for the next refill. */
        const size_t consumed = (bytecnt / GAUSS_RAND_BYTES) * GAUSS_RAND_BYTES;
        const size_t rem = bytecnt - consumed;
        if (rem > 0)
            memmove(buf, buf + consumed, rem);
        bytecnt = rem;
    }
}
