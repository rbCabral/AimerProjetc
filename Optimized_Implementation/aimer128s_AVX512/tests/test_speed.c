// SPDX-License-Identifier: MIT

#include <stddef.h>
#include <stdint.h>
#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include "../api.h"
#include "../common/cpucycles.h"
#include "../common/speed_print.h"
#include "../common/rng.h"
#include "../bench.h"


#define BENCH 50

int main()
{
  
  size_t mlen = 59;
  size_t smlen = 0;
  uint8_t message[mlen];
  uint8_t m[mlen];
  uint8_t pk[CRYPTO_PUBLICKEYBYTES];
  uint8_t sk[CRYPTO_SECRETKEYBYTES];
  uint8_t sm[mlen + CRYPTO_BYTES];

  
  randombytes((unsigned char*)message, 59);

  // BENCH_FUNCTION(crypto_sign_keypair,pk, sk);
  
  // BENCH_FUNCTION(crypto_sign,sm, &smlen, message, mlen, sk);
  BENCH_FUNCTION(crypto_sign_AVX512,sm, &smlen, message, mlen, sk);
  
  
  BENCH_FUNCTION(crypto_sign_open_AVX512,m, &mlen, sm, smlen, pk);
}
