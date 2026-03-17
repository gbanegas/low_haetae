#ifndef API_H
#define API_H
#include "config.h"

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

#define crypto_sign_keypair    HAETAE_NAMESPACE(keypair)
#define crypto_sign_signature  HAETAE_NAMESPACE(signature)
#define crypto_sign_verify     HAETAE_NAMESPACE(verify)


#include "sign.h"
#include "verify.h"
#include "keygen.h"

#endif
