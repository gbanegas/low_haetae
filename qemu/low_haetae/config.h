#ifndef CONFIG_H
#define CONFIG_H

#define HAETAE_MODE 2

/*
* Build-time knobs (override via CFLAGS):
*   -DSAMPLER=1/2        (1 = one-pass hyperball: "fast, but uses large stack arrays", 2 = two-pass low-stack: "lower stack; matches Bos–Renes–Sprenkels style")
*   -DSTREAM_MATRIX=0/1  (1 = streamed A1, no A1[K] on stack)
*   -DFROZEN_A=0/1/2     (0 = full poly tmp, 1 = poly_frozen tmp (512B), 2 = fused sample+multiply, no tmp buffer)                         
*/

#define SAMPLER 1
#define STREAM_MATRIX 1
#define FROZEN_A 2

#define STACK_PROFILE 1
#define STACK_PROFILE_USE_HAL 1

/* Signing memory mode.
 * 1: low-stack signing (stream/reuse; no cs1/cs2/z1/z2 polyvec temporaries)
 * 0: baseline signing
 */
#define SIGN_LOWSTACK 1

#if HAETAE_MODE == 2
#define CRYPTO_ALGNAME "HAETAE2"
#define HAETAE_NAMESPACETOP haetae2
#define HAETAE_NAMESPACE(s) cryptolab_haetae2_##s
#elif HAETAE_MODE == 3
#define CRYPTO_ALGNAME "HAETAE3"
#define HAETAE_NAMESPACETOP haetae3
#define HAETAE_NAMESPACE(s) cryptolab_haetae3_##s
#elif HAETAE_MODE == 5
#define CRYPTO_ALGNAME "HAETAE5"
#define HAETAE_NAMESPACETOP haetae5
#define HAETAE_NAMESPACE(s) cryptolab_haetae5_##s
#endif
#endif


