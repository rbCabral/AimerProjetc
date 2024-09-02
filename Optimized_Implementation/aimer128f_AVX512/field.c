// SPDX-License-Identifier: MIT

#include "field.h"
#include <stddef.h>
#include <stdint.h>
#include <immintrin.h>

#include <stdio.h>

inline __m512i transpose8_512(__m512i x){
  const __m512i C1N =  _mm512_set1_epi64(0xAA55AA55AA55AA55);
  const __m512i C2N =  _mm512_set1_epi64(0x00AA00AA00AA00AA);
  const __m512i C3N =  _mm512_set1_epi64(0xCCCC3333CCCC3333);
  const __m512i C4N =  _mm512_set1_epi64(0x0000CCCC0000CCCC);
  const __m512i C5N =  _mm512_set1_epi64(0xF0F0F0F00F0F0F0F);
  const __m512i C6N =  _mm512_set1_epi64(0x00000000F0F0F0F0);
    // x = _mm512_xor_si512(_mm512_xor_si512(_mm512_and_si512(x,C1N),_mm512_slli_epi64(_mm512_and_si512(x,C2N),7)) ,_mm512_and_si512(_mm512_srli_epi64(x,7) ,C2N));
    // x = _mm512_xor_si512(_mm512_xor_si512(_mm512_and_si512(x,C3N),_mm512_slli_epi64(_mm512_and_si512(x,C4N),14)),_mm512_and_si512(_mm512_srli_epi64(x,14),C4N));
    // x = _mm512_xor_si512(_mm512_xor_si512(_mm512_and_si512(x,C5N),_mm512_slli_epi64(_mm512_and_si512(x,C6N),28)),_mm512_and_si512(_mm512_srli_epi64(x,28),C6N));
    x = _mm512_ternarylogic_epi64(_mm512_and_si512(x,C1N),_mm512_slli_epi64(_mm512_and_si512(x,C2N),7),_mm512_and_si512(_mm512_srli_epi64(x,7) ,C2N),0x96);
    x = _mm512_ternarylogic_epi64(_mm512_and_si512(x,C3N),_mm512_slli_epi64(_mm512_and_si512(x,C4N),14),_mm512_and_si512(_mm512_srli_epi64(x,14) ,C4N),0x96);
    x = _mm512_ternarylogic_epi64(_mm512_and_si512(x,C5N),_mm512_slli_epi64(_mm512_and_si512(x,C6N),28),_mm512_and_si512(_mm512_srli_epi64(x,28) ,C6N),0x96);
    
    return x;
  }

void print_512(__m512i x){
    uint64_t c[8];
    _mm512_storeu_si512((__m512i*)c+0,x);
    printf("%8.16lx %8.16lx %8.16lx %8.16lx %8.16lx %8.16lx %8.16lx %8.16lx\n", c[7],c[6],c[5],c[4],c[3],c[2],c[1],c[0]);
}

void print_256(__m256i x){
    uint64_t c[4];
    _mm256_storeu_si256((__m256i*)c+0,x);
    printf("%8.16lx %8.16lx %8.16lx %8.16lx\n", c[3],c[2],c[1],c[0]);
}

void print_128(__m128i x){
    uint64_t c[2];
    _mm_storeu_si128((__m128i*)c+0,x);
    printf("%8.16lx %8.16lx\n",c[1],c[0]);
}

// static void printBytes(const uint8_t* in, int len){
//     int i;
//     for(i=0;i<len;i++){
//         printf("%.2X", in[i]);
//     }
//     printf("\n");
// }

void GF_to_bytes(uint8_t *out, const GF in)
{
  int i, j;
  for (i = 0; i < AIM2_NUM_WORDS_FIELD; i++)
  {
    uint64_t u = in[i];
    for (j = 0; j < 8; j++)
    {
      *out++ = u;
      u >>= 8;
    }
  }
}

void GF_from_bytes(GF out, const uint8_t *in)
{
  int i, j;
  for (i = 0; i < AIM2_NUM_WORDS_FIELD; i++)
  {
    uint64_t u = 0;
    for (j = 7; j >= 0; j--)
    {
      u = (u << 8) | in[8 * i + j];
    }
    out[i] = u;
  }
}

void GF_set0(GF a)
{
  a[0] = 0;
  a[1] = 0;
}

void GF_copy(GF out, const GF in)
{
  out[0] = in[0];
  out[1] = in[1];
}

void GF_add(GF c, const GF a, const GF b)
{
  c[0] = a[0] ^ b[0];
  c[1] = a[1] ^ b[1];
}

void GF_mul(GF c, const GF a, const GF b)
{
  __m128i temp_a, temp_b;
  __m128i temp_c[2];
  __m128i temp[3];
  __m128i irr = _mm_set_epi64x(0x0, 0x87);

  temp_a = _mm_loadu_si128((const __m128i *)a);
  temp_b = _mm_loadu_si128((const __m128i *)b);

  // polynomial multiplication
  temp_c[0] = _mm_clmulepi64_si128(temp_a, temp_b, 0x00);
  temp_c[1] = _mm_clmulepi64_si128(temp_a, temp_b, 0x11);

  temp[0] = _mm_clmulepi64_si128(temp_a, temp_b, 0x01);
  temp[1] = _mm_clmulepi64_si128(temp_a, temp_b, 0x10);

  temp[0] = _mm_xor_si128(temp[0], temp[1]);
  temp_c[0] = _mm_xor_si128(temp_c[0], _mm_slli_si128(temp[0], 8));
  temp_c[1] = _mm_xor_si128(temp_c[1], _mm_srli_si128(temp[0], 8));

  // modular reduction
  temp[0] = _mm_clmulepi64_si128(temp_c[1], irr, 0x01);
  temp[1] = _mm_slli_si128(temp[0], 8);
  temp[2] = _mm_srli_si128(temp[0], 8);
  temp[2] = _mm_xor_si128(temp[2], temp_c[1]);

  temp_c[0] = _mm_xor_si128(temp_c[0], _mm_clmulepi64_si128(temp[2], irr, 0x00));
  temp_c[0] = _mm_xor_si128(temp_c[0], temp[1]);

  _mm_storeu_si128((__m128i *)c, temp_c[0]);
}

void GF_mul_N_2(GF c[AIMER_N], const GF a[AIMER_N], const GF b)
{
  __m128i x[2], y, z[2], zhi[2], t[4];

  __m128i irr = _mm_set_epi64x(0x0, 0x87);
  y = _mm_loadu_si128((const __m128i *)b);

  for (size_t party = 0; party < AIMER_N; party += 2)
  {
    x[0] = _mm_loadu_si128((const __m128i *)a[party]);
    x[1] = _mm_loadu_si128((const __m128i *)a[party + 1]);

    // polynomial multiplication x2
    t[0] = _mm_clmulepi64_si128(x[0], y, 0x01);
    t[1] = _mm_clmulepi64_si128(x[1], y, 0x01);
    t[2] = _mm_clmulepi64_si128(x[0], y, 0x10);
    t[3] = _mm_clmulepi64_si128(x[1], y, 0x10);

    z[0] = _mm_clmulepi64_si128(x[0], y, 0x00);
    z[1] = _mm_clmulepi64_si128(x[1], y, 0x00);
    zhi[0] = _mm_clmulepi64_si128(x[0], y, 0x11);
    zhi[1] = _mm_clmulepi64_si128(x[1], y, 0x11);

    t[0] = _mm_xor_si128(t[0], t[2]);
    t[1] = _mm_xor_si128(t[1], t[3]);

    z[0] = _mm_xor_si128(z[0], _mm_slli_si128(t[0], 8));
    z[1] = _mm_xor_si128(z[1], _mm_slli_si128(t[1], 8));
    zhi[0] = _mm_xor_si128(zhi[0], _mm_srli_si128(t[0], 8));
    zhi[1] = _mm_xor_si128(zhi[1], _mm_srli_si128(t[1], 8));

    // modular reduction x2
    t[0] = _mm_clmulepi64_si128(zhi[0], irr, 0x01);
    t[1] = _mm_clmulepi64_si128(zhi[1], irr, 0x01);

    z[0] = _mm_xor_si128(z[0], _mm_slli_si128(t[0], 8));
    z[1] = _mm_xor_si128(z[1], _mm_slli_si128(t[1], 8));

    zhi[0] = _mm_xor_si128(zhi[0], _mm_srli_si128(t[0], 8));
    zhi[1] = _mm_xor_si128(zhi[1], _mm_srli_si128(t[1], 8));

    t[0] = _mm_clmulepi64_si128(zhi[0], irr, 0x00);
    t[1] = _mm_clmulepi64_si128(zhi[1], irr, 0x00);

    z[0] = _mm_xor_si128(z[0], t[0]);
    z[1] = _mm_xor_si128(z[1], t[1]);

    _mm_storeu_si128((__m128i *)c[party], z[0]);
    _mm_storeu_si128((__m128i *)c[party + 1], z[1]);
  }
}

void GF_mul_N(GF c[AIMER_N], const GF a[AIMER_N], const GF b)
{
  __m512i X[2], Y, Z[2], Zhi[2], T[4];
  __m512i Irr = _mm512_set_epi64(0x0, 0x87,0x0, 0x87,0x0, 0x87,0x0, 0x87);
  
  Y = _mm512_broadcast_i64x2(_mm_loadu_si128((const __m128i *)b));
  

  for (size_t party = 0; party < AIMER_N; party += 8)
  {
    X[0] = _mm512_loadu_si512((const __m512i *)a[party]);
    X[1] = _mm512_loadu_si512((const __m512i *)a[party]+1);

    // polynomial multiplication x2
    T[0] = _mm512_clmulepi64_epi128(X[0], Y, 0x01);
    T[1] = _mm512_clmulepi64_epi128(X[1], Y, 0x01);
    T[2] = _mm512_clmulepi64_epi128(X[0], Y, 0x10);
    T[3] = _mm512_clmulepi64_epi128(X[1], Y, 0x10);

    Z[0]   = _mm512_clmulepi64_epi128(X[0], Y, 0x00);
    Z[1]   = _mm512_clmulepi64_epi128(X[1], Y, 0x00);
    Zhi[0] = _mm512_clmulepi64_epi128(X[0], Y, 0x11);
    Zhi[1] = _mm512_clmulepi64_epi128(X[1], Y, 0x11);

    T[0] = _mm512_xor_si512(T[0], T[2]);
    T[1] = _mm512_xor_si512(T[1], T[3]);

    Z[0]   = _mm512_xor_si512(Z[0]  , _mm512_bslli_epi128(T[0], 8));
    Z[1]   = _mm512_xor_si512(Z[1]  , _mm512_bslli_epi128(T[1], 8));
    Zhi[0] = _mm512_xor_si512(Zhi[0], _mm512_bsrli_epi128(T[0], 8));
    Zhi[1] = _mm512_xor_si512(Zhi[1], _mm512_bsrli_epi128(T[1], 8));

    // modular reduction x2
    T[0] = _mm512_clmulepi64_epi128(Zhi[0], Irr, 0x01);
    T[1] = _mm512_clmulepi64_epi128(Zhi[1], Irr, 0x01);

    Z[0] = _mm512_xor_si512(Z[0], _mm512_bslli_epi128(T[0], 8));
    Z[1] = _mm512_xor_si512(Z[1], _mm512_bslli_epi128(T[1], 8));

    Zhi[0] = _mm512_xor_si512(Zhi[0], _mm512_bsrli_epi128(T[0], 8));
    Zhi[1] = _mm512_xor_si512(Zhi[1], _mm512_bsrli_epi128(T[1], 8));

    T[0] = _mm512_clmulepi64_epi128(Zhi[0], Irr, 0x00);
    T[1] = _mm512_clmulepi64_epi128(Zhi[1], Irr, 0x00);

    Z[0] = _mm512_xor_si512(Z[0], T[0]);
    Z[1] = _mm512_xor_si512(Z[1], T[1]);

    _mm512_storeu_si512((__m512i *)c[party]    , Z[0]);
    _mm512_storeu_si512((__m512i *)c[party] + 1, Z[1]);
  }
}


void GF_mul_N_unr(GF c[AIMER_N], const GF a[AIMER_N], const GF b)
{
  __m512i X[4], Y, Z[4], Zhi[4], T[8];
  __m512i Irr = _mm512_set_epi64(0x0, 0x87,0x0, 0x87,0x0, 0x87,0x0, 0x87);
  
  Y = _mm512_broadcast_i64x2(_mm_loadu_si128((const __m128i *)b));
  

  // for (size_t party = 0; party < AIMER_N; party += 8)
  // {
  X[0] = _mm512_loadu_si512((const __m512i *)a[0]);
  X[1] = _mm512_loadu_si512((const __m512i *)a[0]+1);
  X[2] = _mm512_loadu_si512((const __m512i *)a[0]+2);
  X[3] = _mm512_loadu_si512((const __m512i *)a[0]+3);

    
  T[0] = _mm512_clmulepi64_epi128(X[0], Y, 0x01);
  T[1] = _mm512_clmulepi64_epi128(X[1], Y, 0x01);
  T[2] = _mm512_clmulepi64_epi128(X[0], Y, 0x10);
  T[3] = _mm512_clmulepi64_epi128(X[1], Y, 0x10);
  T[4] = _mm512_clmulepi64_epi128(X[2], Y, 0x01);
  T[5] = _mm512_clmulepi64_epi128(X[3], Y, 0x01);
  T[6] = _mm512_clmulepi64_epi128(X[2], Y, 0x10);
  T[7] = _mm512_clmulepi64_epi128(X[3], Y, 0x10);

  Z[0]   = _mm512_clmulepi64_epi128(X[0], Y, 0x00);
  Z[1]   = _mm512_clmulepi64_epi128(X[1], Y, 0x00);
  Zhi[0] = _mm512_clmulepi64_epi128(X[0], Y, 0x11);
  Zhi[1] = _mm512_clmulepi64_epi128(X[1], Y, 0x11);
  Z[2]   = _mm512_clmulepi64_epi128(X[2], Y, 0x00);
  Z[3]   = _mm512_clmulepi64_epi128(X[3], Y, 0x00);
  Zhi[2] = _mm512_clmulepi64_epi128(X[2], Y, 0x11);
  Zhi[3] = _mm512_clmulepi64_epi128(X[3], Y, 0x11);

  T[0] = _mm512_xor_si512(T[0], T[2]);
  T[1] = _mm512_xor_si512(T[1], T[3]);
  T[2] = _mm512_xor_si512(T[4], T[6]);
  T[3] = _mm512_xor_si512(T[5], T[7]);

  Z[0]   = _mm512_xor_si512(Z[0]  , _mm512_bslli_epi128(T[0], 8));
  Z[1]   = _mm512_xor_si512(Z[1]  , _mm512_bslli_epi128(T[1], 8));
  Zhi[0] = _mm512_xor_si512(Zhi[0], _mm512_bsrli_epi128(T[0], 8));
  Zhi[1] = _mm512_xor_si512(Zhi[1], _mm512_bsrli_epi128(T[1], 8));
  Z[2]   = _mm512_xor_si512(Z[2]  , _mm512_bslli_epi128(T[2], 8));
  Z[3]   = _mm512_xor_si512(Z[3]  , _mm512_bslli_epi128(T[3], 8));
  Zhi[2] = _mm512_xor_si512(Zhi[2], _mm512_bsrli_epi128(T[2], 8));
  Zhi[3] = _mm512_xor_si512(Zhi[3], _mm512_bsrli_epi128(T[3], 8));

    // modular reduction x2
  T[0] = _mm512_clmulepi64_epi128(Zhi[0], Irr, 0x01);
  T[1] = _mm512_clmulepi64_epi128(Zhi[1], Irr, 0x01);
  T[2] = _mm512_clmulepi64_epi128(Zhi[2], Irr, 0x01);
  T[3] = _mm512_clmulepi64_epi128(Zhi[3], Irr, 0x01);

  Z[0] = _mm512_xor_si512(Z[0], _mm512_bslli_epi128(T[0], 8));
  Z[1] = _mm512_xor_si512(Z[1], _mm512_bslli_epi128(T[1], 8));
  Z[2] = _mm512_xor_si512(Z[2], _mm512_bslli_epi128(T[2], 8));
  Z[3] = _mm512_xor_si512(Z[3], _mm512_bslli_epi128(T[3], 8));

  Zhi[0] = _mm512_xor_si512(Zhi[0], _mm512_bsrli_epi128(T[0], 8));
  Zhi[1] = _mm512_xor_si512(Zhi[1], _mm512_bsrli_epi128(T[1], 8));
  Zhi[2] = _mm512_xor_si512(Zhi[2], _mm512_bsrli_epi128(T[2], 8));
  Zhi[3] = _mm512_xor_si512(Zhi[3], _mm512_bsrli_epi128(T[3], 8));

  T[0] = _mm512_clmulepi64_epi128(Zhi[0], Irr, 0x00);
  T[1] = _mm512_clmulepi64_epi128(Zhi[1], Irr, 0x00);
  T[2] = _mm512_clmulepi64_epi128(Zhi[2], Irr, 0x00);
  T[3] = _mm512_clmulepi64_epi128(Zhi[3], Irr, 0x00);

  Z[0] = _mm512_xor_si512(Z[0], T[0]);
  Z[1] = _mm512_xor_si512(Z[1], T[1]);
  Z[2] = _mm512_xor_si512(Z[2], T[2]);
  Z[3] = _mm512_xor_si512(Z[3], T[3]);

  _mm512_storeu_si512((__m512i *)c[0]    , Z[0]);
  _mm512_storeu_si512((__m512i *)c[0] + 1, Z[1]);
  _mm512_storeu_si512((__m512i *)c[0] + 2, Z[2]);
  _mm512_storeu_si512((__m512i *)c[0] + 3, Z[3]);
  
}

void GF_mul_add(GF c, const GF a, const GF b)
{
  __m128i temp_a, temp_b;
  __m128i temp_c[2];
  __m128i temp[3];
  __m128i irr = _mm_set_epi64x(0x0, 0x87);

  temp_a = _mm_loadu_si128((const __m128i *)a);
  temp_b = _mm_loadu_si128((const __m128i *)b);

  // polynomial multiplication
  temp_c[0] = _mm_clmulepi64_si128(temp_a, temp_b, 0x00);
  temp_c[1] = _mm_clmulepi64_si128(temp_a, temp_b, 0x11);

  temp[0] = _mm_clmulepi64_si128(temp_a, temp_b, 0x01);
  temp[1] = _mm_clmulepi64_si128(temp_a, temp_b, 0x10);

  temp[0] = _mm_xor_si128(temp[0], temp[1]);
  temp_c[0] = _mm_xor_si128(temp_c[0], _mm_slli_si128(temp[0], 8));
  temp_c[1] = _mm_xor_si128(temp_c[1], _mm_srli_si128(temp[0], 8));

  // modular reduction
  temp[0] = _mm_clmulepi64_si128(temp_c[1], irr, 0x01);
  temp[1] = _mm_slli_si128(temp[0], 8);
  temp[2] = _mm_srli_si128(temp[0], 8);
  temp[2] = _mm_xor_si128(temp[2], temp_c[1]);

  temp_c[0] = _mm_xor_si128(temp_c[0], _mm_clmulepi64_si128(temp[2], irr, 0x00));
  temp_c[0] = _mm_xor_si128(temp_c[0], temp[1]);

  temp[0] = _mm_loadu_si128((const __m128i *)c);
  temp_c[0] = _mm_xor_si128(temp_c[0], temp[0]);
  _mm_storeu_si128((__m128i *)c, temp_c[0]);
}


void GF_mul_add_N(GF c[AIMER_N], const GF a[AIMER_N], const GF b)
{
  __m512i X[2], Y, Z[2], Zhi[2], T[4];
  __m512i Irr = _mm512_set_epi64(0x0, 0x87,0x0, 0x87,0x0, 0x87,0x0, 0x87);
  
  Y = _mm512_broadcast_i64x2(_mm_loadu_si128((const __m128i *)b));
  

  for (size_t party = 0; party < AIMER_N; party += 8)
  {
    X[0] = _mm512_loadu_si512((const __m512i *)a[party]);
    X[1] = _mm512_loadu_si512((const __m512i *)a[party]+1);

    // polynomial multiplication x2
    T[0] = _mm512_clmulepi64_epi128(X[0], Y, 0x01);
    T[1] = _mm512_clmulepi64_epi128(X[1], Y, 0x01);
    T[2] = _mm512_clmulepi64_epi128(X[0], Y, 0x10);
    T[3] = _mm512_clmulepi64_epi128(X[1], Y, 0x10);

    Z[0]   = _mm512_clmulepi64_epi128(X[0], Y, 0x00);
    Z[1]   = _mm512_clmulepi64_epi128(X[1], Y, 0x00);
    Zhi[0] = _mm512_clmulepi64_epi128(X[0], Y, 0x11);
    Zhi[1] = _mm512_clmulepi64_epi128(X[1], Y, 0x11);

    T[0] = _mm512_xor_si512(T[0], T[2]);
    T[1] = _mm512_xor_si512(T[1], T[3]);

    Z[0]   = _mm512_xor_si512(Z[0]  , _mm512_bslli_epi128(T[0], 8));
    Z[1]   = _mm512_xor_si512(Z[1]  , _mm512_bslli_epi128(T[1], 8));
    Zhi[0] = _mm512_xor_si512(Zhi[0], _mm512_bsrli_epi128(T[0], 8));
    Zhi[1] = _mm512_xor_si512(Zhi[1], _mm512_bsrli_epi128(T[1], 8));

    // load c
    X[0] = _mm512_loadu_si512((const __m512i *)c[party]);
    X[1] = _mm512_loadu_si512((const __m512i *)c[party]+1);

    // modular reduction x2
    T[0] = _mm512_clmulepi64_epi128(Zhi[0], Irr, 0x01);
    T[1] = _mm512_clmulepi64_epi128(Zhi[1], Irr, 0x01);

    Z[0] = _mm512_xor_si512(Z[0], _mm512_bslli_epi128(T[0], 8));
    Z[1] = _mm512_xor_si512(Z[1], _mm512_bslli_epi128(T[1], 8));
    Zhi[0] = _mm512_xor_si512(Zhi[0], _mm512_bsrli_epi128(T[0], 8));
    Zhi[1] = _mm512_xor_si512(Zhi[1], _mm512_bsrli_epi128(T[1], 8));

    T[0] = _mm512_clmulepi64_epi128(Zhi[0], Irr, 0x00);
    T[1] = _mm512_clmulepi64_epi128(Zhi[1], Irr, 0x00);

    Z[0] = _mm512_xor_si512(Z[0], X[0]);
    Z[1] = _mm512_xor_si512(Z[1], X[1]);

    Z[0] = _mm512_xor_si512(Z[0], T[0]);
    Z[1] = _mm512_xor_si512(Z[1], T[1]);

    _mm512_storeu_si512((__m512i *)c[party]    , Z[0]);
    _mm512_storeu_si512((__m512i *)c[party] + 1, Z[1]);
  }
}

void GF_mul_add_N_2(GF c[AIMER_N], const GF a[AIMER_N], const GF b)
{
  __m128i x[2], y, z[2], zhi[2], t[4];
  __m128i irr = _mm_set_epi64x(0x0, 0x87);
  y = _mm_loadu_si128((const __m128i *)b);

  for (size_t party = 0; party < AIMER_N; party += 2)
  {
    x[0] = _mm_loadu_si128((const __m128i *)a[party]);
    x[1] = _mm_loadu_si128((const __m128i *)a[party + 1]);

    // polynomial multiplication x2
    t[0] = _mm_clmulepi64_si128(x[0], y, 0x01);
    t[1] = _mm_clmulepi64_si128(x[1], y, 0x01);
    t[2] = _mm_clmulepi64_si128(x[0], y, 0x10);
    t[3] = _mm_clmulepi64_si128(x[1], y, 0x10);

    z[0] = _mm_clmulepi64_si128(x[0], y, 0x00);
    z[1] = _mm_clmulepi64_si128(x[1], y, 0x00);
    zhi[0] = _mm_clmulepi64_si128(x[0], y, 0x11);
    zhi[1] = _mm_clmulepi64_si128(x[1], y, 0x11);

    t[0] = _mm_xor_si128(t[0], t[2]);
    t[1] = _mm_xor_si128(t[1], t[3]);

    z[0] = _mm_xor_si128(z[0], _mm_slli_si128(t[0], 8));
    z[1] = _mm_xor_si128(z[1], _mm_slli_si128(t[1], 8));
    zhi[0] = _mm_xor_si128(zhi[0], _mm_srli_si128(t[0], 8));
    zhi[1] = _mm_xor_si128(zhi[1], _mm_srli_si128(t[1], 8));

    // load c
    x[0] = _mm_loadu_si128((const __m128i *)c[party]);
    x[1] = _mm_loadu_si128((const __m128i *)c[party + 1]);

    // modular reduction x2
    t[0] = _mm_clmulepi64_si128(zhi[0], irr, 0x01);
    t[1] = _mm_clmulepi64_si128(zhi[1], irr, 0x01);

    z[0] = _mm_xor_si128(z[0], _mm_slli_si128(t[0], 8));
    z[1] = _mm_xor_si128(z[1], _mm_slli_si128(t[1], 8));
    zhi[0] = _mm_xor_si128(zhi[0], _mm_srli_si128(t[0], 8));
    zhi[1] = _mm_xor_si128(zhi[1], _mm_srli_si128(t[1], 8));

    t[0] = _mm_clmulepi64_si128(zhi[0], irr, 0x00);
    t[1] = _mm_clmulepi64_si128(zhi[1], irr, 0x00);

    z[0] = _mm_xor_si128(z[0], x[0]);
    z[1] = _mm_xor_si128(z[1], x[1]);

    z[0] = _mm_xor_si128(z[0], t[0]);
    z[1] = _mm_xor_si128(z[1], t[1]);

    _mm_storeu_si128((__m128i *)c[party], z[0]);
    _mm_storeu_si128((__m128i *)c[party + 1], z[1]);
  }
}

void GF_sqr(GF c, const GF a)
{
  __m128i temp_a;
  __m128i temp_c[2];
  __m128i temp[3];
  __m128i irr = _mm_set_epi64x(0x0, 0x87);

  temp_a = _mm_loadu_si128((const __m128i *)a);

  // polynomial squaring
  temp_c[0] = _mm_clmulepi64_si128(temp_a, temp_a, 0x00);
  temp_c[1] = _mm_clmulepi64_si128(temp_a, temp_a, 0x11);

  // modular reduction
  temp[0] = _mm_clmulepi64_si128(temp_c[1], irr, 0x01);
  temp[1] = _mm_slli_si128(temp[0], 8);
  temp[2] = _mm_srli_si128(temp[0], 8);
  temp[2] = _mm_xor_si128(temp[2], temp_c[1]);

  temp_c[0] = _mm_xor_si128(temp_c[0], _mm_clmulepi64_si128(temp[2], irr, 0x00));
  temp_c[0] = _mm_xor_si128(temp_c[0], temp[1]);

  _mm_storeu_si128((__m128i *)c, temp_c[0]);
}


void GF_sqr_N(GF c[AIMER_N], const GF a[AIMER_N])
{
  __m512i Z[2], Zhi[2], T[2];
  __m512i Irr = _mm512_set_epi64(0x0, 0x87,0x0, 0x87,0x0, 0x87,0x0, 0x87);

  for (size_t party = 0; party < AIMER_N; party += 8)
  {
    T[0] = _mm512_loadu_si512((const __m512i *)a[party]);
    T[1] = _mm512_loadu_si512((const __m512i *)a[party]+1);

    // polynomial squaring x4
    Z[0]   = _mm512_clmulepi64_epi128(T[0], T[0], 0x00);
    Z[1]   = _mm512_clmulepi64_epi128(T[1], T[1], 0x00);

    Zhi[0] = _mm512_clmulepi64_epi128(T[0], T[0], 0x11);
    Zhi[1] = _mm512_clmulepi64_epi128(T[1], T[1], 0x11);

    // modular reduction x4
    T[0] = _mm512_clmulepi64_epi128(Zhi[0], Irr, 0x01);
    T[1] = _mm512_clmulepi64_epi128(Zhi[1], Irr, 0x01);

    Z[0] = _mm512_xor_si512(Z[0], _mm512_bslli_epi128(T[0], 8));
    Z[1] = _mm512_xor_si512(Z[1], _mm512_bslli_epi128(T[1], 8));
    Zhi[0] = _mm512_xor_si512(Zhi[0], _mm512_bsrli_epi128(T[0], 8));
    Zhi[1] = _mm512_xor_si512(Zhi[1], _mm512_bsrli_epi128(T[1], 8));

    T[0] = _mm512_clmulepi64_epi128(Zhi[0], Irr, 0x00);
    T[1] = _mm512_clmulepi64_epi128(Zhi[1], Irr, 0x00);

    Z[0] = _mm512_xor_si512(Z[0], T[0]);
    Z[1] = _mm512_xor_si512(Z[1], T[1]);

    _mm512_storeu_si512((__m512i *)c[party]    , Z[0]);
    _mm512_storeu_si512((__m512i *)c[party] + 1, Z[1]);
  }
}

void GF_sqr_N_256(GF c[AIMER_N], const GF a[AIMER_N])
{
  __m128i z[4], zhi[4], t[4];
  __m128i irr = _mm_set_epi64x(0x0, 0x87);

  for (size_t party = 0; party < AIMER_N; party += 4)
  {
    t[0] = _mm_loadu_si128((const __m128i *)a[party]);
    t[1] = _mm_loadu_si128((const __m128i *)a[party + 1]);
    t[2] = _mm_loadu_si128((const __m128i *)a[party + 2]);
    t[3] = _mm_loadu_si128((const __m128i *)a[party + 3]);

    // polynomial squaring x4
    z[0] = _mm_clmulepi64_si128(t[0], t[0], 0x00);
    z[1] = _mm_clmulepi64_si128(t[1], t[1], 0x00);
    z[2] = _mm_clmulepi64_si128(t[2], t[2], 0x00);
    z[3] = _mm_clmulepi64_si128(t[3], t[3], 0x00);

    zhi[0] = _mm_clmulepi64_si128(t[0], t[0], 0x11);
    zhi[1] = _mm_clmulepi64_si128(t[1], t[1], 0x11);
    zhi[2] = _mm_clmulepi64_si128(t[2], t[2], 0x11);
    zhi[3] = _mm_clmulepi64_si128(t[3], t[3], 0x11);

    // modular reduction x4
    t[0] = _mm_clmulepi64_si128(zhi[0], irr, 0x01);
    t[1] = _mm_clmulepi64_si128(zhi[1], irr, 0x01);
    t[2] = _mm_clmulepi64_si128(zhi[2], irr, 0x01);
    t[3] = _mm_clmulepi64_si128(zhi[3], irr, 0x01);

    z[0] = _mm_xor_si128(z[0], _mm_slli_si128(t[0], 8));
    z[1] = _mm_xor_si128(z[1], _mm_slli_si128(t[1], 8));
    zhi[0] = _mm_xor_si128(zhi[0], _mm_srli_si128(t[0], 8));
    zhi[1] = _mm_xor_si128(zhi[1], _mm_srli_si128(t[1], 8));

    z[2] = _mm_xor_si128(z[2], _mm_slli_si128(t[2], 8));
    z[3] = _mm_xor_si128(z[3], _mm_slli_si128(t[3], 8));
    zhi[2] = _mm_xor_si128(zhi[2], _mm_srli_si128(t[2], 8));
    zhi[3] = _mm_xor_si128(zhi[3], _mm_srli_si128(t[3], 8));

    t[0] = _mm_clmulepi64_si128(zhi[0], irr, 0x00);
    t[1] = _mm_clmulepi64_si128(zhi[1], irr, 0x00);
    t[2] = _mm_clmulepi64_si128(zhi[2], irr, 0x00);
    t[3] = _mm_clmulepi64_si128(zhi[3], irr, 0x00);

    z[0] = _mm_xor_si128(z[0], t[0]);
    z[1] = _mm_xor_si128(z[1], t[1]);
    z[2] = _mm_xor_si128(z[2], t[2]);
    z[3] = _mm_xor_si128(z[3], t[3]);

    _mm_storeu_si128((__m128i *)c[party], z[0]);
    _mm_storeu_si128((__m128i *)c[party + 1], z[1]);
    _mm_storeu_si128((__m128i *)c[party + 2], z[2]);
    _mm_storeu_si128((__m128i *)c[party + 3], z[3]);
  }
}




void GF_transposed_matmul(GF c, const GF a, const GF b[AIM2_NUM_BITS_FIELD])
{
  unsigned int i, j;
  const uint32_t *a_ptr = (uint32_t *)a;
  const __m256i shift = _mm256_set_epi32(0, 2, 4, 6, 1, 3, 5, 7);
  __m256i temp[2] = {_mm256_setzero_si256(), _mm256_setzero_si256()};
  __m256i matrix_data[4];


  for (i = 0; i < 4; i++)
  {
    __m256i index = _mm256_set1_epi32(a_ptr[i]);
    for (j = 32 * (i + 1); j > 32 * i; j -= 8)
    {
      __m256i mask = _mm256_sllv_epi32(index, shift);
      matrix_data[0] = _mm256_loadu_si256((const __m256i *)&b[j - 2]);
      matrix_data[1] = _mm256_loadu_si256((const __m256i *)&b[j - 4]);
      matrix_data[2] = _mm256_loadu_si256((const __m256i *)&b[j - 6]);
      matrix_data[3] = _mm256_loadu_si256((const __m256i *)&b[j - 8]);

      temp[0] = _mm256_xor_si256(temp[0],
        _mm256_and_si256(matrix_data[0],
          _mm256_srai_epi32(_mm256_shuffle_epi32(mask, 0xff), 31)));
      temp[1] = _mm256_xor_si256(temp[1],
        _mm256_and_si256(matrix_data[1],
          _mm256_srai_epi32(_mm256_shuffle_epi32(mask, 0xaa), 31)));
      temp[0] = _mm256_xor_si256(temp[0],
        _mm256_and_si256(matrix_data[2],
          _mm256_srai_epi32(_mm256_shuffle_epi32(mask, 0x55), 31)));
      temp[1] = _mm256_xor_si256(temp[1],
        _mm256_and_si256(matrix_data[3],
          _mm256_srai_epi32(_mm256_shuffle_epi32(mask, 0x00), 31)));
      index = _mm256_slli_epi32(index, 8);
    }
  }

  temp[0] = _mm256_xor_si256(temp[0], temp[1]);

  _mm_storeu_si128((__m128i *)c,
                  _mm_xor_si128(_mm256_extracti128_si256(temp[0], 0),
                  _mm256_extracti128_si256(temp[0], 1)));
}

// void GF_transposed_matmul_2(GF c[AIM2_NUM_BITS_FIELD], const GF a[AIM2_NUM_BITS_FIELD], const GF b[AIM2_NUM_BITS_FIELD])
// {
//   uint64_t matrix_U_l[128];
//   uint64_t matrix_U_h[128];
//   int k=0;
//   /*posso mudar na geração para já ficar nessa organização*/
//     for(int i=0;i<16;i++){
//       for(int j=i*8;j<i*8+8;j++){
//         matrix_U_l[j]   = a[k][0];
//         matrix_U_h[j]   = a[k][1];
//         k++;
//       }
//     }
  
//     const __m256i zero = _mm256_setzero_si256();
//   __m256i m0, m1, m2, m3, t;
//   const __m512i one = _mm512_set1_epi64(0x01);

//   uint8_t stateOut0[16];  
//   uint8_t stateOut1[16];  
//   uint8_t o;
 
//     for(int j=0,id=0;j<128;j+=2,id++){
      
//       __m512i Al0,Al1,Ah0,Ah1,Tmp, Tmp2, aux;

//       Tmp = _mm512_loadu_si512((const __m512i *)b[j]);
      
//       Al0 =  _mm512_broadcast_i32x2(_mm512_castsi512_si128(Tmp));
//       Ah0 =  _mm512_broadcast_i32x2(_mm512_castsi512_si128(_mm512_bsrli_epi128(Tmp,8)));
      
//       Tmp = _mm512_loadu_si512((const __m512i *)b[j+1]);
//       Al1 =  _mm512_broadcast_i32x2(_mm512_castsi512_si128(Tmp));
//       Ah1 =  _mm512_broadcast_i32x2(_mm512_castsi512_si128(_mm512_bsrli_epi128(Tmp,8)));
      
//       for(int i=0;i<8;i++){
//         Tmp  = _mm512_loadu_si512((const __m512i *)matrix_U_l+i);
//         Tmp2 = _mm512_loadu_si512((const __m512i *)matrix_U_h+i);

//         aux = _mm512_xor_epi64(_mm512_and_epi64(Al0,Tmp),_mm512_and_epi64(Ah0,Tmp2));
//         aux = _mm512_popcnt_epi64(aux);
//         aux = _mm512_and_epi64(aux,one);
//         o = _mm512_cmp_epi64_mask (aux, one, 0x00);
//         stateOut0[i] = o;

//         aux = _mm512_xor_epi64(_mm512_and_epi64(Al1,Tmp),_mm512_and_epi64(Ah1,Tmp2));
//         aux = _mm512_popcnt_epi64(aux);
//         aux = _mm512_and_epi64(aux,one);
//         o = _mm512_cmp_epi64_mask (aux, one, 0x00);
//         stateOut0[i+8] = o;
//       }
      
//       for(int i=8,k1=0;i<16;i++,k1++){
//         Tmp  = _mm512_loadu_si512((const __m512i *)matrix_U_l+i);
//         Tmp2 = _mm512_loadu_si512((const __m512i *)matrix_U_h+i);

//         aux = _mm512_xor_epi64(_mm512_and_epi64(Al0,Tmp),_mm512_and_epi64(Ah0,Tmp2));
//         aux = _mm512_popcnt_epi64(aux);
//         aux = _mm512_and_epi64(aux,one);
//         o = _mm512_cmp_epi64_mask (aux, one, 0x00);
//         stateOut1[k1] = o;

//         aux = _mm512_xor_epi64(_mm512_and_epi64(Al1,Tmp),_mm512_and_epi64(Ah1,Tmp2));
//         aux = _mm512_popcnt_epi64(aux);
//         aux = _mm512_and_epi64(aux,one);
//         o = _mm512_cmp_epi64_mask (aux, one, 0x00);
//         stateOut1[k1+8] = o;
//       }
      
//       GF_from_bytes(c[id], stateOut0);
//       GF_from_bytes(c[id+64], stateOut1);
//     }
// }


void GF_transposed_matmul_add_N_AA(GF c[AIMER_N], const GF a[AIMER_N],const GF b[AIM2_NUM_BITS_FIELD])
{
  ALIGN2 uint8_t tmpL[64];
  ALIGN2 uint8_t tmpH[64];
  ALIGN2 uint8_t tmpL1[64];
  ALIGN2 uint8_t tmpH1[64];
  
  __m512i Al,Ah, Aux0,Aux1, TmpL, TmpH, a0,a1,a01,a11,Aux01,Aux11;
 
  const __m512i m0 =  _mm512_set_epi8(63,55,47,39,31,23,15,7,
                                      62,54,46,38,30,22,14,6,
                                      61,53,45,37,29,21,13,5,
                                      60,52,44,36,28,20,12,4,
                                      59,51,43,35,27,19,11,3,
                                      58,50,42,34,26,18,10,2,
                                      57,49,41,33,25,17,9,1,
                                      56,48,40,32,24,16,8,0);

   
  Aux0 = _mm512_loadu_si512((const __m512i *)a);
  Aux1 = _mm512_loadu_si512((const __m512i *)a+1);

  TmpL = _mm512_unpacklo_epi64 (Aux0,Aux1);
  TmpH = _mm512_unpackhi_epi64 (Aux0,Aux1);

  TmpL = _mm512_permutexvar_epi8(m0,TmpL);
  TmpH = _mm512_permutexvar_epi8(m0,TmpH);

  Al = _mm512_loadu_si512((const __m512i *)c);
  Ah = _mm512_loadu_si512((const __m512i *)c+1);
    
  TmpL = transpose8_512(TmpL);
  TmpH = transpose8_512(TmpH);
  _mm512_store_si512((__m512i*)tmpL,TmpL);
  _mm512_store_si512((__m512i*)tmpH,TmpH);

  Aux0 = _mm512_unpacklo_epi64 (Al,Ah);
  Aux1 = _mm512_unpackhi_epi64 (Al,Ah);

  //------------------
  Aux01 = _mm512_loadu_si512((const __m512i *)a+2);
  Aux11 = _mm512_loadu_si512((const __m512i *)a+3);

  TmpL = _mm512_unpacklo_epi64 (Aux01,Aux11);
  TmpH = _mm512_unpackhi_epi64 (Aux01,Aux11);

  TmpL = _mm512_permutexvar_epi8(m0,TmpL);
  TmpH = _mm512_permutexvar_epi8(m0,TmpH);

  Al = _mm512_loadu_si512((const __m512i *)c+2);
  Ah = _mm512_loadu_si512((const __m512i *)c+3);
    
  TmpL = transpose8_512(TmpL);
  TmpH = transpose8_512(TmpH);
  _mm512_store_si512((__m512i*)tmpL1,TmpL);
  _mm512_store_si512((__m512i*)tmpH1,TmpH);

  Aux01 = _mm512_unpacklo_epi64 (Al,Ah);
  Aux11 = _mm512_unpackhi_epi64 (Al,Ah);
  
  //------------------
  

  for(int j=0,k=64;j<64;j++,k++){
    a0 = _mm512_set1_epi64(b[j][0]);
    a1 = _mm512_set1_epi64(b[j][1]);
    a01 = _mm512_set1_epi64(b[k][0]);
    a11 = _mm512_set1_epi64(b[k][1]);

    Aux0 = _mm512_mask_xor_epi64(Aux0,tmpL[j],a0,Aux0);
    Aux0 = _mm512_mask_xor_epi64(Aux0,tmpH[j],a01,Aux0);
    Aux1 = _mm512_mask_xor_epi64(Aux1,tmpL[j],a1,Aux1);
    Aux1 = _mm512_mask_xor_epi64(Aux1,tmpH[j],a11,Aux1);

    Aux01 = _mm512_mask_xor_epi64(Aux01,tmpL1[j],a0,Aux01);
    Aux01 = _mm512_mask_xor_epi64(Aux01,tmpH1[j],a01,Aux01);
    Aux11 = _mm512_mask_xor_epi64(Aux11,tmpL1[j],a1,Aux11);
    Aux11 = _mm512_mask_xor_epi64(Aux11,tmpH1[j],a11,Aux11);
  }

  Al = _mm512_unpacklo_epi64 (Aux0,Aux1);
  Ah = _mm512_unpackhi_epi64 (Aux0,Aux1);
  Aux0 = _mm512_unpacklo_epi64 (Aux01,Aux11);
  Aux1 = _mm512_unpackhi_epi64 (Aux01,Aux11);

  _mm512_storeu_si512((__m512i*)c  ,Al);
  _mm512_storeu_si512((__m512i*)c+1,Ah);
  _mm512_storeu_si512((__m512i*)c+2,Aux0);
  _mm512_storeu_si512((__m512i*)c+3,Aux1);

}

void GF_transposed_matmul_add_N_A(GF c[AIMER_N], const GF a[AIMER_N],const uint64_t bl[AIM2_NUM_BITS_FIELD],const uint64_t bh[AIM2_NUM_BITS_FIELD])
{
  ALIGN2 uint8_t tmpL[64];
  ALIGN2 uint8_t tmpH[64];
  ALIGN2 uint8_t tmpL1[64];
  ALIGN2 uint8_t tmpH1[64];
  
  __m512i Al,Ah, Aux0,Aux1, TmpL, TmpH, a0,a1,a01,a11,Aux01,Aux11;
 
  const __m512i m0 =  _mm512_set_epi8(63,55,47,39,31,23,15,7,
                                      62,54,46,38,30,22,14,6,
                                      61,53,45,37,29,21,13,5,
                                      60,52,44,36,28,20,12,4,
                                      59,51,43,35,27,19,11,3,
                                      58,50,42,34,26,18,10,2,
                                      57,49,41,33,25,17,9,1,
                                      56,48,40,32,24,16,8,0);

   
  Aux0 = _mm512_loadu_si512((const __m512i *)a);
  Aux1 = _mm512_loadu_si512((const __m512i *)a+1);

  TmpL = _mm512_unpacklo_epi64 (Aux0,Aux1);
  TmpH = _mm512_unpackhi_epi64 (Aux0,Aux1);

  TmpL = _mm512_permutexvar_epi8(m0,TmpL);
  TmpH = _mm512_permutexvar_epi8(m0,TmpH);

  Al = _mm512_loadu_si512((const __m512i *)c);
  Ah = _mm512_loadu_si512((const __m512i *)c+1);
    
  TmpL = transpose8_512(TmpL);
  TmpH = transpose8_512(TmpH);
  _mm512_store_si512((__m512i*)tmpL,TmpL);
  _mm512_store_si512((__m512i*)tmpH,TmpH);

  Aux0 = _mm512_unpacklo_epi64 (Al,Ah);
  Aux1 = _mm512_unpackhi_epi64 (Al,Ah);

  //------------------
  Aux01 = _mm512_loadu_si512((const __m512i *)a+2);
  Aux11 = _mm512_loadu_si512((const __m512i *)a+3);

  TmpL = _mm512_unpacklo_epi64 (Aux01,Aux11);
  TmpH = _mm512_unpackhi_epi64 (Aux01,Aux11);

  TmpL = _mm512_permutexvar_epi8(m0,TmpL);
  TmpH = _mm512_permutexvar_epi8(m0,TmpH);

  Al = _mm512_loadu_si512((const __m512i *)c+2);
  Ah = _mm512_loadu_si512((const __m512i *)c+3);
    
  TmpL = transpose8_512(TmpL);
  TmpH = transpose8_512(TmpH);
  _mm512_store_si512((__m512i*)tmpL1,TmpL);
  _mm512_store_si512((__m512i*)tmpH1,TmpH);

  Aux01 = _mm512_unpacklo_epi64 (Al,Ah);
  Aux11 = _mm512_unpackhi_epi64 (Al,Ah);
  
  //------------------
  

  for(int j=0,k=64;j<64;j++,k++){
    a0 = _mm512_set1_epi64(bl[j]);
    a1 = _mm512_set1_epi64(bh[j]);
    a01 = _mm512_set1_epi64(bl[k]);
    a11 = _mm512_set1_epi64(bh[k]);

    Aux0 = _mm512_mask_xor_epi64(Aux0,tmpL[j],a0,Aux0);
    Aux0 = _mm512_mask_xor_epi64(Aux0,tmpH[j],a01,Aux0);
    Aux1 = _mm512_mask_xor_epi64(Aux1,tmpL[j],a1,Aux1);
    Aux1 = _mm512_mask_xor_epi64(Aux1,tmpH[j],a11,Aux1);

    Aux01 = _mm512_mask_xor_epi64(Aux01,tmpL1[j],a0,Aux01);
    Aux01 = _mm512_mask_xor_epi64(Aux01,tmpH1[j],a01,Aux01);
    Aux11 = _mm512_mask_xor_epi64(Aux11,tmpL1[j],a1,Aux11);
    Aux11 = _mm512_mask_xor_epi64(Aux11,tmpH1[j],a11,Aux11);
  }

  Al = _mm512_unpacklo_epi64 (Aux0,Aux1);
  Ah = _mm512_unpackhi_epi64 (Aux0,Aux1);
  Aux0 = _mm512_unpacklo_epi64 (Aux01,Aux11);
  Aux1 = _mm512_unpackhi_epi64 (Aux01,Aux11);

  _mm512_storeu_si512((__m512i*)c  ,Al);
  _mm512_storeu_si512((__m512i*)c+1,Ah);
  _mm512_storeu_si512((__m512i*)c+2,Aux0);
  _mm512_storeu_si512((__m512i*)c+3,Aux1);

}


void GF_transposed_matmul_add_N_A_v12(GF c[AIMER_N], const GF a[AIMER_N],const uint64_t bl[AIM2_NUM_BITS_FIELD],const uint64_t bh[AIM2_NUM_BITS_FIELD])
{
  ALIGN2 uint8_t tmpL[64];
  ALIGN2 uint8_t tmpH[64];
  ALIGN2 uint8_t tmpL1[64];
  ALIGN2 uint8_t tmpH1[64];
  
  __m512i Al,Ah, Aux0,Aux1, TmpL, TmpH, a0,a1,a01,a11,Aux01,Aux11;
 
  const __m512i m0 =  _mm512_set_epi8(63,55,47,39,31,23,15,7,
                                      62,54,46,38,30,22,14,6,
                                      61,53,45,37,29,21,13,5,
                                      60,52,44,36,28,20,12,4,
                                      59,51,43,35,27,19,11,3,
                                      58,50,42,34,26,18,10,2,
                                      57,49,41,33,25,17,9,1,
                                      56,48,40,32,24,16,8,0);

   
  Aux0 = _mm512_loadu_si512((const __m512i *)a);
  Aux1 = _mm512_loadu_si512((const __m512i *)a+1);

  TmpL = _mm512_unpacklo_epi64 (Aux0,Aux1);
  TmpH = _mm512_unpackhi_epi64 (Aux0,Aux1);

  TmpL = _mm512_permutexvar_epi8(m0,TmpL);
  TmpH = _mm512_permutexvar_epi8(m0,TmpH);

  Al = _mm512_loadu_si512((const __m512i *)c);
  Ah = _mm512_loadu_si512((const __m512i *)c+1);
    
  TmpL = transpose8_512(TmpL);
  TmpH = transpose8_512(TmpH);
  _mm512_store_si512((__m512i*)tmpL,TmpL);
  _mm512_store_si512((__m512i*)tmpH,TmpH);

  Aux0 = _mm512_unpacklo_epi64 (Al,Ah);
  Aux1 = _mm512_unpackhi_epi64 (Al,Ah);

  //------------------
  Aux01 = _mm512_loadu_si512((const __m512i *)a+2);
  Aux11 = _mm512_loadu_si512((const __m512i *)a+3);

  TmpL = _mm512_unpacklo_epi64 (Aux01,Aux11);
  TmpH = _mm512_unpackhi_epi64 (Aux01,Aux11);

  TmpL = _mm512_permutexvar_epi8(m0,TmpL);
  TmpH = _mm512_permutexvar_epi8(m0,TmpH);

  Al = _mm512_loadu_si512((const __m512i *)c+2);
  Ah = _mm512_loadu_si512((const __m512i *)c+3);
    
  TmpL = transpose8_512(TmpL);
  TmpH = transpose8_512(TmpH);
  _mm512_store_si512((__m512i*)tmpL1,TmpL);
  _mm512_store_si512((__m512i*)tmpH1,TmpH);

  Aux01 = _mm512_unpacklo_epi64 (Al,Ah);
  Aux11 = _mm512_unpackhi_epi64 (Al,Ah);
  
  //------------------
  

  for(int j=0,k=64;j<64;j++,k++){
    a0 = _mm512_set1_epi64(bl[j]);
    a1 = _mm512_set1_epi64(bh[j]);
    a01 = _mm512_set1_epi64(bl[k]);
    a11 = _mm512_set1_epi64(bh[k]);

    Aux0 = _mm512_mask_xor_epi64(Aux0,tmpL[j],a0,Aux0);
    Aux0 = _mm512_mask_xor_epi64(Aux0,tmpH[j],a01,Aux0);
    Aux1 = _mm512_mask_xor_epi64(Aux1,tmpL[j],a1,Aux1);
    Aux1 = _mm512_mask_xor_epi64(Aux1,tmpH[j],a11,Aux1);

    Aux01 = _mm512_mask_xor_epi64(Aux01,tmpL1[j],a0,Aux01);
    Aux01 = _mm512_mask_xor_epi64(Aux01,tmpH1[j],a01,Aux01);
    Aux11 = _mm512_mask_xor_epi64(Aux11,tmpL1[j],a1,Aux11);
    Aux11 = _mm512_mask_xor_epi64(Aux11,tmpH1[j],a11,Aux11);
  }

  Al = _mm512_unpacklo_epi64 (Aux0,Aux1);
  Ah = _mm512_unpackhi_epi64 (Aux0,Aux1);
  Aux0 = _mm512_unpacklo_epi64 (Aux01,Aux11);
  Aux1 = _mm512_unpackhi_epi64 (Aux01,Aux11);

  _mm512_storeu_si512((__m512i*)c  ,Al);
  _mm512_storeu_si512((__m512i*)c+1,Ah);
  _mm512_storeu_si512((__m512i*)c+2,Aux0);
  _mm512_storeu_si512((__m512i*)c+3,Aux1);

}

void GF_transposed_matmul_add_N_A_v11(GF c[AIMER_N], const GF a[AIMER_N],const uint64_t bl[AIM2_NUM_BITS_FIELD],const uint64_t bh[AIM2_NUM_BITS_FIELD])
{
  ALIGN2 uint8_t tmpL[64];
  ALIGN2 uint8_t tmpH[64];
  
  __m512i Al,Ah, Aux0,Aux1, TmpL, TmpH, a0,a1,a01,a11;
 
  const __m512i m0 =  _mm512_set_epi8(63,55,47,39,31,23,15,7,
                                      62,54,46,38,30,22,14,6,
                                      61,53,45,37,29,21,13,5,
                                      60,52,44,36,28,20,12,4,
                                      59,51,43,35,27,19,11,3,
                                      58,50,42,34,26,18,10,2,
                                      57,49,41,33,25,17,9,1,
                                      56,48,40,32,24,16,8,0);

  
  for(int p=0;p<4;p+=2){
    Aux0 = _mm512_loadu_si512((const __m512i *)a+p);
    Aux1 = _mm512_loadu_si512((const __m512i *)a+p+1);

    TmpL = _mm512_unpacklo_epi64 (Aux0,Aux1);
    TmpH = _mm512_unpackhi_epi64 (Aux0,Aux1);

    TmpL = _mm512_permutexvar_epi8(m0,TmpL);
    TmpH = _mm512_permutexvar_epi8(m0,TmpH);

    Al = _mm512_loadu_si512((const __m512i *)c+p);
    Ah = _mm512_loadu_si512((const __m512i *)c+p+1);
    
    TmpL = transpose8_512(TmpL);
    TmpH = transpose8_512(TmpH);
    _mm512_store_si512((__m512i*)tmpL,TmpL);
    _mm512_store_si512((__m512i*)tmpH,TmpH);

    Aux0 = _mm512_unpacklo_epi64 (Al,Ah);
    Aux1 = _mm512_unpackhi_epi64 (Al,Ah);


    for(int j=0,k=64;j<64;j++,k++){
      a0 = _mm512_set1_epi64(bl[j]);
      a1 = _mm512_set1_epi64(bh[j]);
      a01 = _mm512_set1_epi64(bl[k]);
      a11 = _mm512_set1_epi64(bh[k]);

      Aux0 = _mm512_mask_xor_epi64(Aux0,tmpL[j],a0,Aux0);
      Aux0 = _mm512_mask_xor_epi64(Aux0,tmpH[j],a01,Aux0);
      Aux1 = _mm512_mask_xor_epi64(Aux1,tmpL[j],a1,Aux1);
      Aux1 = _mm512_mask_xor_epi64(Aux1,tmpH[j],a11,Aux1);
    }

    Al = _mm512_unpacklo_epi64 (Aux0,Aux1);
    Ah = _mm512_unpackhi_epi64 (Aux0,Aux1);

    _mm512_storeu_si512((__m512i*)c+p  ,Al);
    _mm512_storeu_si512((__m512i*)c+p+1,Ah);
  }  

}


void GF_transposed_matmul_add_N_A_v10(GF c[AIMER_N], const GF a[AIMER_N],const uint64_t bl[AIM2_NUM_BITS_FIELD],const uint64_t bh[AIM2_NUM_BITS_FIELD])
{

  
  __m256i Al,Ah,Al1,Ah1, Aux0,Aux1,Aux2,Aux3, TmpL, TmpH,TmpL1, TmpH1,r0,r1,r2,r3;
  const __m256i one =  _mm256_set1_epi64x(1);
  const __m256i zero = _mm256_setzero_si256();
  
  for(int p=0;p<8;p+=4){
    Aux0 = _mm256_loadu_si256((const __m256i *)a+p);
    Aux1 = _mm256_loadu_si256((const __m256i *)a+p+1);
    Aux2 = _mm256_loadu_si256((const __m256i *)a+p+2);
    Aux3 = _mm256_loadu_si256((const __m256i *)a+p+3);

    TmpL  = _mm256_unpacklo_epi64 (Aux0,Aux1);
    TmpH  = _mm256_unpackhi_epi64 (Aux0,Aux1);
    TmpL1 = _mm256_unpacklo_epi64 (Aux2,Aux3);
    TmpH1 = _mm256_unpackhi_epi64 (Aux2,Aux3);

    Al  = _mm256_loadu_si256((const __m256i *)c+p);
    Ah  = _mm256_loadu_si256((const __m256i *)c+p+1);
    Al1 = _mm256_loadu_si256((const __m256i *)c+p+2);
    Ah1 = _mm256_loadu_si256((const __m256i *)c+p+3);

    r0 = _mm256_unpacklo_epi64(Al,Ah);
    r1 = _mm256_unpackhi_epi64(Al,Ah);
    r2 = _mm256_unpacklo_epi64(Al1,Ah1);
    r3 = _mm256_unpackhi_epi64(Al1,Ah1);

    for(int j=0,k=64;j<64;j++,k++){
      Aux0 = _mm256_set1_epi64x(bl[j]);
      Aux1 = _mm256_set1_epi64x(bh[j]);
      Aux2 = _mm256_set1_epi64x(bl[k]);
      Aux3 = _mm256_set1_epi64x(bh[k]);   
      
      Al  = _mm256_and_si256(TmpL,one);
      Ah  = _mm256_and_si256(TmpH,one);
      Al1 = _mm256_and_si256(TmpL1,one);
      Ah1 = _mm256_and_si256(TmpH1,one);

      Al  = _mm256_sub_epi64(zero,Al);
      Ah  = _mm256_sub_epi64(zero,Ah);
      Al1 = _mm256_sub_epi64(zero,Al1);
      Ah1 = _mm256_sub_epi64(zero,Ah1);

      r0 = _mm256_xor_si256(_mm256_and_si256(Al,Aux0),r0);
      r1 = _mm256_xor_si256(_mm256_and_si256(Al,Aux1),r1);
      r0 = _mm256_xor_si256(_mm256_and_si256(Ah,Aux2),r0);
      r1 = _mm256_xor_si256(_mm256_and_si256(Ah,Aux3),r1);

      r2 = _mm256_xor_si256(_mm256_and_si256(Al1,Aux0),r2);
      r3 = _mm256_xor_si256(_mm256_and_si256(Al1,Aux1),r3);
      r2 = _mm256_xor_si256(_mm256_and_si256(Ah1,Aux2),r2);
      r3 = _mm256_xor_si256(_mm256_and_si256(Ah1,Aux3),r3);

      TmpL  = _mm256_srli_epi64(TmpL,1);
      TmpH  = _mm256_srli_epi64(TmpH,1);
      TmpL1 = _mm256_srli_epi64(TmpL1,1);
      TmpH1 = _mm256_srli_epi64(TmpH1,1);
    }
    
    Al  = _mm256_unpacklo_epi64 (r0,r1);
    Ah  = _mm256_unpackhi_epi64 (r0,r1);
    Al1 = _mm256_unpacklo_epi64 (r2,r3);
    Ah1 = _mm256_unpackhi_epi64 (r2,r3);
    
    _mm256_storeu_si256((__m256i*)c+p  ,Al);
    _mm256_storeu_si256((__m256i*)c+p+1,Ah);
    _mm256_storeu_si256((__m256i*)c+p+2,Al1);
    _mm256_storeu_si256((__m256i*)c+p+3,Ah1);
  }  
}

void GF_transposed_matmul_add_N_A_v9(GF c[AIMER_N], const GF a[AIMER_N],const uint64_t bl[AIM2_NUM_BITS_FIELD],const uint64_t bh[AIM2_NUM_BITS_FIELD])
{

 
  __m512i Al,Ah,Al1,Ah1, Aux0,Aux1, TmpL, TmpH, TmpL1, TmpH1, r0,r1,r2,r3;
  const __m512i one =  _mm512_set1_epi64(1);
  const __m512i zero = _mm512_setzero_si512();
  
    Al = _mm512_loadu_si512((const __m512i *)a);
    Ah = _mm512_loadu_si512((const __m512i *)a+1);
    TmpL = _mm512_unpacklo_epi64 (Al,Ah);
    TmpH = _mm512_unpackhi_epi64 (Al,Ah);

    Al1 = _mm512_loadu_si512((const __m512i *)a+2);
    Ah1 = _mm512_loadu_si512((const __m512i *)a+3);
    TmpL1 = _mm512_unpacklo_epi64 (Al1,Ah1);
    TmpH1 = _mm512_unpackhi_epi64 (Al1,Ah1);

    Al = _mm512_loadu_si512((const __m512i *)c);
    Ah = _mm512_loadu_si512((const __m512i *)c+1);
    r0 = _mm512_unpacklo_epi64 (Al,Ah);
    r1 = _mm512_unpackhi_epi64 (Al,Ah);

    Al1 = _mm512_loadu_si512((const __m512i *)c+2);
    Ah1 = _mm512_loadu_si512((const __m512i *)c+3);
    r2 = _mm512_unpacklo_epi64 (Al1,Ah1);
    r3 = _mm512_unpackhi_epi64 (Al1,Ah1);

    
    for(int j=0,k=64;j<64;j++,k++){
      Aux0 = _mm512_set1_epi64(bl[j]);
      Aux1 = _mm512_set1_epi64(bh[j]);

      Al = _mm512_and_si512(TmpL,one);
      Al1 = _mm512_and_si512(TmpL1,one);
      Al = _mm512_sub_epi64(zero,Al);
      Al1 = _mm512_sub_epi64(zero,Al1);

      r0 = _mm512_ternarylogic_epi64(Al,Aux0,r0,0x6a);
      r1 = _mm512_ternarylogic_epi64(Al,Aux1,r1,0x6a);
      r2 = _mm512_ternarylogic_epi64(Al1,Aux0,r2,0x6a);
      r3 = _mm512_ternarylogic_epi64(Al1,Aux1,r3,0x6a);
      // r0 = _mm512_xor_si512(_mm512_and_si512(Al,Aux0),r0);
      // r1 = _mm512_xor_si512(_mm512_and_si512(Al,Aux1),r1);
      // r2 = _mm512_xor_si512(_mm512_and_si512(Al1,Aux0),r2);
      // r3 = _mm512_xor_si512(_mm512_and_si512(Al1,Aux1),r3);
      TmpL = _mm512_srli_epi64(TmpL,1);
      TmpL1 = _mm512_srli_epi64(TmpL1,1);

      Aux0 = _mm512_set1_epi64(bl[k]);
      Aux1 = _mm512_set1_epi64(bh[k]);
      
      Al = _mm512_and_si512(TmpH,one);
      Al = _mm512_sub_epi64(zero,Al);
      Al1 = _mm512_and_si512(TmpH1,one);
      Al1 = _mm512_sub_epi64(zero,Al1);

      r0 = _mm512_ternarylogic_epi64(Al,Aux0,r0,0x6a);
      r1 = _mm512_ternarylogic_epi64(Al,Aux1,r1,0x6a);
      r2 = _mm512_ternarylogic_epi64(Al1,Aux0,r2,0x6a);
      r3 = _mm512_ternarylogic_epi64(Al1,Aux1,r3,0x6a);
      // r0 = _mm512_xor_si512(_mm512_and_si512(Al,Aux0),r0);
      // r1 = _mm512_xor_si512(_mm512_and_si512(Al,Aux1),r1);
      // r2 = _mm512_xor_si512(_mm512_and_si512(Al1,Aux0),r2);
      // r3 = _mm512_xor_si512(_mm512_and_si512(Al1,Aux1),r3);
      TmpH = _mm512_srli_epi64(TmpH,1);
      TmpH1 = _mm512_srli_epi64(TmpH1,1);
    }

    Al  = _mm512_unpacklo_epi64 (r0,r1);
    Ah  = _mm512_unpackhi_epi64 (r0,r1);
    Al1 = _mm512_unpacklo_epi64 (r2,r3);
    Ah1 = _mm512_unpackhi_epi64 (r2,r3);
    
    _mm512_storeu_si512((__m512i*)c  ,Al);
    _mm512_storeu_si512((__m512i*)c+1,Ah);
    _mm512_storeu_si512((__m512i*)c+2,Al1);
    _mm512_storeu_si512((__m512i*)c+3,Ah1);
   
}


void GF_transposed_matmul_add_N_A_v8(GF c[AIMER_N], const GF a[AIMER_N],const uint64_t bl[AIM2_NUM_BITS_FIELD],const uint64_t bh[AIM2_NUM_BITS_FIELD])
{
  ALIGN2 uint8_t tmpL[64];
  ALIGN2 uint8_t tmpH[64];
  
  __m512i Al,Ah, Aux0,Aux1, TmpL, TmpH;
 
  const __m512i m0 =  _mm512_set_epi8(63,55,47,39,31,23,15,7,
                                      62,54,46,38,30,22,14,6,
                                      61,53,45,37,29,21,13,5,
                                      60,52,44,36,28,20,12,4,
                                      59,51,43,35,27,19,11,3,
                                      58,50,42,34,26,18,10,2,
                                      57,49,41,33,25,17,9,1,
                                      56,48,40,32,24,16,8,0);

  
  for(int p=0;p<4;p+=2){
    Aux0 = _mm512_loadu_si512((const __m512i *)a+p);
    Aux1 = _mm512_loadu_si512((const __m512i *)a+p+1);

    TmpL = _mm512_unpacklo_epi64 (Aux0,Aux1);
    TmpH = _mm512_unpackhi_epi64 (Aux0,Aux1);

    TmpL = _mm512_permutexvar_epi8(m0,TmpL);
    TmpH = _mm512_permutexvar_epi8(m0,TmpH);

    Al = _mm512_loadu_si512((const __m512i *)c+p);
    Ah = _mm512_loadu_si512((const __m512i *)c+p+1);
    
    TmpL = transpose8_512(TmpL);
    TmpH = transpose8_512(TmpH);
    _mm512_store_si512((__m512i*)tmpL,TmpL);
    _mm512_store_si512((__m512i*)tmpH,TmpH);

    Aux0 = _mm512_unpacklo_epi64 (Al,Ah);
    Aux1 = _mm512_unpackhi_epi64 (Al,Ah);

    for(int j=0,k=64;j<64;j++,k++){
      Aux0 = _mm512_ternarylogic_epi64(Aux0,_mm512_maskz_set1_epi64(tmpL[j],bl[j]),_mm512_maskz_set1_epi64(tmpH[j],bl[k]),0x96);
      Aux1 = _mm512_ternarylogic_epi64(Aux1,_mm512_maskz_set1_epi64(tmpL[j],bh[j]),_mm512_maskz_set1_epi64(tmpH[j],bh[k]),0x96);
    }

    Al = _mm512_unpacklo_epi64 (Aux0,Aux1);
    Ah = _mm512_unpackhi_epi64 (Aux0,Aux1);

    _mm512_storeu_si512((__m512i*)c+p  ,Al);
    _mm512_storeu_si512((__m512i*)c+p+1,Ah);
  }  

}

void GF_transposed_matmul_add_N_A_v7(GF c[AIMER_N], const GF a[AIMER_N],const uint64_t bl[AIM2_NUM_BITS_FIELD],const uint64_t bh[AIM2_NUM_BITS_FIELD])
{

  __m512i Al,Ah,Al1,Ah1, Aux0,Aux1, TmpL, TmpH, TmpL1, TmpH1, r0,r1,r2,r3;
  
    Al = _mm512_loadu_si512((const __m512i *)a);
    Ah = _mm512_loadu_si512((const __m512i *)a+1);
    Al1 = _mm512_loadu_si512((const __m512i *)a+2);
    Ah1 = _mm512_loadu_si512((const __m512i *)a+3);

    TmpL = _mm512_unpacklo_epi64 (Al,Ah);
    TmpH = _mm512_unpackhi_epi64 (Al,Ah);
    TmpL1 = _mm512_unpacklo_epi64 (Al1,Ah1);
    TmpH1 = _mm512_unpackhi_epi64 (Al1,Ah1);

    const __m512i one =  _mm512_set1_epi64(1);
    const __m512i zero = _mm512_setzero_si512();
    
    Al = _mm512_loadu_si512((const __m512i *)c);
    Ah = _mm512_loadu_si512((const __m512i *)c+1);
    Al1 = _mm512_loadu_si512((const __m512i *)c+2);
    Ah1 = _mm512_loadu_si512((const __m512i *)c+3);

    r0 = _mm512_unpacklo_epi64 (Al,Ah);
    r1 = _mm512_unpackhi_epi64 (Al,Ah);
    r2 = _mm512_unpacklo_epi64 (Al1,Ah1);
    r3 = _mm512_unpackhi_epi64 (Al1,Ah1);

    
    for(int j=0;j<64;j++){
      Aux0 = _mm512_set1_epi64(bl[j]);
      Aux1 = _mm512_set1_epi64(bh[j]);
      
      Al = _mm512_and_si512(TmpL,one);
      Al = _mm512_sub_epi64(zero,Al);
      Al1 = _mm512_and_si512(TmpL1,one);
      Al1 = _mm512_sub_epi64(zero,Al1);

      r0 = _mm512_xor_si512(_mm512_and_si512(Al,Aux0),r0);
      r1 = _mm512_xor_si512(_mm512_and_si512(Al,Aux1),r1);
      r2 = _mm512_xor_si512(_mm512_and_si512(Al1,Aux0),r2);
      r3 = _mm512_xor_si512(_mm512_and_si512(Al1,Aux1),r3);
      TmpL = _mm512_srli_epi64(TmpL,1);
      TmpL1 = _mm512_srli_epi64(TmpL1,1);
    }

    for(int j=64;j<128;j++){
      Aux0 = _mm512_set1_epi64(bl[j]);
      Aux1 = _mm512_set1_epi64(bh[j]);
      
      Al = _mm512_and_si512(TmpH,one);
      Al = _mm512_sub_epi64(zero,Al);
      Al1 = _mm512_and_si512(TmpH1,one);
      Al1 = _mm512_sub_epi64(zero,Al1);

      r0 = _mm512_xor_si512(_mm512_and_si512(Al,Aux0),r0);
      r1 = _mm512_xor_si512(_mm512_and_si512(Al,Aux1),r1);
      r2 = _mm512_xor_si512(_mm512_and_si512(Al1,Aux0),r2);
      r3 = _mm512_xor_si512(_mm512_and_si512(Al1,Aux1),r3);
      TmpH = _mm512_srli_epi64(TmpH,1);
      TmpH1 = _mm512_srli_epi64(TmpH1,1);
    }

    Al  = _mm512_unpacklo_epi64 (r0,r1);
    Ah  = _mm512_unpackhi_epi64 (r0,r1);
    Al1 = _mm512_unpacklo_epi64 (r2,r3);
    Ah1 = _mm512_unpackhi_epi64 (r2,r3);
    
    _mm512_storeu_si512((__m512i*)c  ,Al);
    _mm512_storeu_si512((__m512i*)c+1,Ah);
    _mm512_storeu_si512((__m512i*)c+2,Al1);
    _mm512_storeu_si512((__m512i*)c+3,Ah1);
   
}


void GF_transposed_matmul_add_N_A_v6(GF c[AIMER_N], const GF a[AIMER_N],const uint64_t bl[AIM2_NUM_BITS_FIELD],const uint64_t bh[AIM2_NUM_BITS_FIELD])
{
 
  __m512i Al,Ah, Aux0,Aux1, TmpL, TmpH,r0,r1;
  
  for(int p=0;p<4;p+=2){
    Aux0 = _mm512_loadu_si512((const __m512i *)a+p);
    Aux1 = _mm512_loadu_si512((const __m512i *)a+p+1);

    TmpL = _mm512_unpacklo_epi64 (Aux0,Aux1);
    TmpH = _mm512_unpackhi_epi64 (Aux0,Aux1);

    const __m512i one =  _mm512_set1_epi64(1);
    const __m512i zero = _mm512_setzero_si512();
    // r0 = _mm512_setzero_si512();//_mm512_loadu_si512((const __m512i *)c+p);
    // r1 = _mm512_setzero_si512();//_mm512_loadu_si512((const __m512i *)c+p+1);
    Al = _mm512_loadu_si512((const __m512i *)c+p);
    Ah = _mm512_loadu_si512((const __m512i *)c+p+1);

    r0 = _mm512_unpacklo_epi64 (Al,Ah);
    r1 = _mm512_unpackhi_epi64 (Al,Ah);

    
    for(int j=0;j<64;j++){
      Aux0 = _mm512_set1_epi64(bl[j]);
      Aux1 = _mm512_set1_epi64(bh[j]);
      
      Al = _mm512_and_si512(TmpL,one);
      Al = _mm512_sub_epi64(zero,Al);

      r0 = _mm512_xor_si512(_mm512_and_si512(Al,Aux0),r0);
      r1 = _mm512_xor_si512(_mm512_and_si512(Al,Aux1),r1);
      TmpL = _mm512_srli_epi64(TmpL,1);
    }

    for(int j=64;j<128;j++){
      Aux0 = _mm512_set1_epi64(bl[j]);
      Aux1 = _mm512_set1_epi64(bh[j]);
      
      Al = _mm512_and_si512(TmpH,one);
      Al = _mm512_sub_epi64(zero,Al);

      r0 = _mm512_xor_si512(_mm512_and_si512(Al,Aux0),r0);
      r1 = _mm512_xor_si512(_mm512_and_si512(Al,Aux1),r1);
      TmpH = _mm512_srli_epi64(TmpH,1);
    }
    Al = _mm512_unpacklo_epi64 (r0,r1);
    Ah = _mm512_unpackhi_epi64 (r0,r1);
    
    _mm512_storeu_si512((__m512i*)c+p  ,Al);
    _mm512_storeu_si512((__m512i*)c+p+1,Ah);
  }  
}


void GF_transposed_matmul_add_N_A_v5(GF c[AIMER_N], const GF a[AIMER_N],const uint64_t bl[AIM2_NUM_BITS_FIELD],const uint64_t bh[AIM2_NUM_BITS_FIELD])
{
  uint8_t tmpL[64];
  uint8_t tmpH[64];
  uint8_t tmpL1[64];
  uint8_t tmpH1[64];
  
  __m512i Al,Ah, Aux0,Aux1, TmpL, TmpH;
  __m512i Al1,Ah1, Aux2,Aux3, TmpL1, TmpH1;
  
  const __m512i m0 =  _mm512_set_epi8(63,55,47,39,31,23,15,7,
                                      62,54,46,38,30,22,14,6,
                                      61,53,45,37,29,21,13,5,
                                      60,52,44,36,28,20,12,4,
                                      59,51,43,35,27,19,11,3,
                                      58,50,42,34,26,18,10,2,
                                      57,49,41,33,25,17,9,1,
                                      56,48,40,32,24,16,8,0);

  
  
  Aux0 = _mm512_loadu_si512((const __m512i *)a+0);
  Aux1 = _mm512_loadu_si512((const __m512i *)a+1);
  Aux2 = _mm512_loadu_si512((const __m512i *)a+2);
  Aux3 = _mm512_loadu_si512((const __m512i *)a+3);

  TmpL  = _mm512_unpacklo_epi64 (Aux0,Aux1);
  TmpH  = _mm512_unpackhi_epi64 (Aux0,Aux1);
  TmpL1 = _mm512_unpacklo_epi64 (Aux2,Aux3);
  TmpH1 = _mm512_unpackhi_epi64 (Aux2,Aux3);

  TmpL  = _mm512_permutexvar_epi8(m0,TmpL);
  TmpH  = _mm512_permutexvar_epi8(m0,TmpH);
  TmpL1 = _mm512_permutexvar_epi8(m0,TmpL1);
  TmpH1 = _mm512_permutexvar_epi8(m0,TmpH1);

  TmpL = transpose8_512(TmpL);
  TmpH = transpose8_512(TmpH);
  TmpL1 = transpose8_512(TmpL1);
  TmpH1 = transpose8_512(TmpH1);
  _mm512_storeu_si512((__m512i*)tmpL,TmpL);
  _mm512_storeu_si512((__m512i*)tmpH,TmpH);
  _mm512_storeu_si512((__m512i*)tmpL1,TmpL1);
  _mm512_storeu_si512((__m512i*)tmpH1,TmpH1);

  Aux0 = _mm512_maskz_set1_epi64(tmpL[0],bl[0]);
  Aux1 = _mm512_maskz_set1_epi64(tmpL[0],bh[0]);    
  Aux2 = _mm512_maskz_set1_epi64(tmpL1[0],bl[0]);
  Aux3 = _mm512_maskz_set1_epi64(tmpL1[0],bh[0]);    
  for(int j=1;j<64;j++){
    Aux0 = _mm512_xor_epi64(Aux0,_mm512_maskz_set1_epi64(tmpL[j],bl[j]));
    Aux1 = _mm512_xor_epi64(Aux1,_mm512_maskz_set1_epi64(tmpL[j],bh[j]));
    Aux2 = _mm512_xor_epi64(Aux2,_mm512_maskz_set1_epi64(tmpL1[j],bl[j]));
    Aux3 = _mm512_xor_epi64(Aux3,_mm512_maskz_set1_epi64(tmpL1[j],bh[j]));
    }
    
  for(int j=64,k=0;j<128;j++,k++){
    Aux0 = _mm512_xor_epi64(Aux0,_mm512_maskz_set1_epi64(tmpH[k],bl[j]));
    Aux1 = _mm512_xor_epi64(Aux1,_mm512_maskz_set1_epi64(tmpH[k],bh[j]));
    Aux2 = _mm512_xor_epi64(Aux2,_mm512_maskz_set1_epi64(tmpH1[k],bl[j]));
    Aux3 = _mm512_xor_epi64(Aux3,_mm512_maskz_set1_epi64(tmpH1[k],bh[j]));
  }

  TmpL  = _mm512_loadu_si512((const __m512i *)c+0);
  TmpH  = _mm512_loadu_si512((const __m512i *)c+1);
  TmpL1 = _mm512_loadu_si512((const __m512i *)c+2);
  TmpH1 = _mm512_loadu_si512((const __m512i *)c+3);

  Al = _mm512_unpacklo_epi64 (TmpL,TmpH);
  Ah = _mm512_unpackhi_epi64 (TmpL,TmpH);
  Al1 = _mm512_unpacklo_epi64 (TmpL1,TmpH1);
  Ah1 = _mm512_unpackhi_epi64 (TmpL1,TmpH1);

  Aux0 = _mm512_xor_epi64(Aux0, Al);
  Aux1 = _mm512_xor_epi64(Aux1, Ah);
  Aux2 = _mm512_xor_epi64(Aux2, Al1);
  Aux3 = _mm512_xor_epi64(Aux3, Ah1);

  Al  = _mm512_unpacklo_epi64 (Aux0,Aux1);
  Ah  = _mm512_unpackhi_epi64 (Aux0,Aux1);
  Al1 = _mm512_unpacklo_epi64 (Aux2,Aux3);
  Ah1 = _mm512_unpackhi_epi64 (Aux2,Aux3);

  _mm512_storeu_si512((__m512i*)c+0,Al);
  _mm512_storeu_si512((__m512i*)c+1,Ah);
  _mm512_storeu_si512((__m512i*)c+2,Al1);
  _mm512_storeu_si512((__m512i*)c+3,Ah1);

}

void GF_transposed_matmul_add_N_A_v4(GF c[AIMER_N], const GF a[AIMER_N],const uint64_t bl[AIM2_NUM_BITS_FIELD],const uint64_t bh[AIM2_NUM_BITS_FIELD])
{
  

  uint8_t tmpL[64];
  uint8_t tmpH[64];
  
  __m512i Al,Ah, Aux0,Aux1, TmpL, TmpH;
  
  
  const __m512i m0 =  _mm512_set_epi8(63,55,47,39,31,23,15,7,
                                      62,54,46,38,30,22,14,6,
                                      61,53,45,37,29,21,13,5,
                                      60,52,44,36,28,20,12,4,
                                      59,51,43,35,27,19,11,3,
                                      58,50,42,34,26,18,10,2,
                                      57,49,41,33,25,17,9,1,
                                      56,48,40,32,24,16,8,0);

  
  for(int p=0;p<4;p+=2){
    Aux0 = _mm512_loadu_si512((const __m512i *)a+p);
    Aux1 = _mm512_loadu_si512((const __m512i *)a+p+1);

    TmpL = _mm512_unpacklo_epi64 (Aux0,Aux1);
    TmpH = _mm512_unpackhi_epi64 (Aux0,Aux1);

    TmpL = _mm512_permutexvar_epi8(m0,TmpL);
    TmpH = _mm512_permutexvar_epi8(m0,TmpH);

    TmpL = transpose8_512(TmpL);
    TmpH = transpose8_512(TmpH);
    _mm512_storeu_si512((__m512i*)tmpL,TmpL);
    _mm512_storeu_si512((__m512i*)tmpH,TmpH);

    Aux0 = _mm512_maskz_set1_epi64(tmpL[0],bl[0]);
    Aux1 = _mm512_maskz_set1_epi64(tmpL[0],bh[0]);
    for(int j=1;j<64;j++){
      Aux0 = _mm512_xor_epi64(Aux0,_mm512_maskz_set1_epi64(tmpL[j],bl[j]));
      Aux1 = _mm512_xor_epi64(Aux1,_mm512_maskz_set1_epi64(tmpL[j],bh[j]));
    }
    
    for(int j=64,k=0;j<128;j++,k++){
      Aux0 = _mm512_xor_epi64(Aux0,_mm512_maskz_set1_epi64(tmpH[k],bl[j]));
      Aux1 = _mm512_xor_epi64(Aux1,_mm512_maskz_set1_epi64(tmpH[k],bh[j]));
    }

    TmpL = _mm512_loadu_si512((const __m512i *)c+p);
    TmpH = _mm512_loadu_si512((const __m512i *)c+p+1);

    Al = _mm512_unpacklo_epi64 (TmpL,TmpH);
    Ah = _mm512_unpackhi_epi64 (TmpL,TmpH);

    Aux0 = _mm512_xor_epi64(Aux0, Al);
    Aux1 = _mm512_xor_epi64(Aux1, Ah);

    Al = _mm512_unpacklo_epi64 (Aux0,Aux1);
    Ah = _mm512_unpackhi_epi64 (Aux0,Aux1);

    
    _mm512_storeu_si512((__m512i*)c+p  ,Al);
    _mm512_storeu_si512((__m512i*)c+p+1,Ah);
  }  

}
  
void GF_transposed_matmul_add_N_A_v3(GF c[AIMER_N], const GF a[AIMER_N],const uint64_t bl[AIM2_NUM_BITS_FIELD],const uint64_t bh[AIM2_NUM_BITS_FIELD])
{
  const __m512i one = _mm512_set1_epi64(0x01);

  uint8_t stateOut0[16];
  uint8_t stateOut1[16];
  uint8_t stateOut2[16];
  uint8_t stateOut3[16];
  
  __m512i Al[4],Ah[4], Aux, Tmp, Tmp2,aux[4];
  __m128i aux2;
  
  for(int j=0;j<16;j+=4){
    Aux = _mm512_loadu_si512((const __m512i *)a[j]);
    aux2 = _mm512_extracti32x4_epi32(Aux,0);
    Al[0] =  _mm512_broadcast_i32x2(aux2);
    Ah[0] =  _mm512_broadcast_i32x2(_mm_bsrli_si128(aux2,8));
    aux2 = _mm512_extracti32x4_epi32(Aux,1);
    Al[1] =  _mm512_broadcast_i32x2(aux2);
    Ah[1] =  _mm512_broadcast_i32x2(_mm_bsrli_si128(aux2,8));
    aux2 = _mm512_extracti32x4_epi32(Aux,2);
    Al[2] =  _mm512_broadcast_i32x2(aux2);
    Ah[2] =  _mm512_broadcast_i32x2(_mm_bsrli_si128(aux2,8));
    aux2 = _mm512_extracti32x4_epi32(Aux,3);
    Al[3] =  _mm512_broadcast_i32x2(aux2);
    Ah[3] =  _mm512_broadcast_i32x2(_mm_bsrli_si128(aux2,8));
    
    for(int i=15,k=0;i>=0;i--,k++){
      Tmp  = _mm512_loadu_si512((const __m512i *)bl+i);
      Tmp2 = _mm512_loadu_si512((const __m512i *)bh+i);
      
      aux[0] = _mm512_xor_epi64(_mm512_and_epi64(Al[0],Tmp),_mm512_and_epi64(Ah[0],Tmp2));
      aux[1] = _mm512_xor_epi64(_mm512_and_epi64(Al[1],Tmp),_mm512_and_epi64(Ah[1],Tmp2));
      aux[2] = _mm512_xor_epi64(_mm512_and_epi64(Al[2],Tmp),_mm512_and_epi64(Ah[2],Tmp2));
      aux[3] = _mm512_xor_epi64(_mm512_and_epi64(Al[3],Tmp),_mm512_and_epi64(Ah[3],Tmp2));
      aux[0] = _mm512_popcnt_epi64(aux[0]);
      aux[1] = _mm512_popcnt_epi64(aux[1]);
      aux[2] = _mm512_popcnt_epi64(aux[2]);
      aux[3] = _mm512_popcnt_epi64(aux[3]);
      aux[0] = _mm512_and_epi64(aux[0],one);
      aux[1] = _mm512_and_epi64(aux[1],one);
      aux[2] = _mm512_and_epi64(aux[2],one);
      aux[3] = _mm512_and_epi64(aux[3],one);
      uint8_t o0 = _mm512_cmp_epi64_mask (aux[0], one, 0x00);
      uint8_t o1 = _mm512_cmp_epi64_mask (aux[1], one, 0x00);
      uint8_t o2 = _mm512_cmp_epi64_mask (aux[2], one, 0x00);
      uint8_t o3 = _mm512_cmp_epi64_mask (aux[3], one, 0x00);
      stateOut0[k] = o0;
      stateOut1[k] = o1;
      stateOut2[k] = o2;
      stateOut3[k] = o3;
    }
    GF r0,r1,r2,r3;
    GF_from_bytes(r0,stateOut0);
    GF_from_bytes(r1,stateOut1);
    GF_from_bytes(r2,stateOut2);
    GF_from_bytes(r3,stateOut3);
    GF_add(c[j+0],c[j+0],r0);
    GF_add(c[j+1],c[j+1],r1);
    GF_add(c[j+2],c[j+2],r2);
    GF_add(c[j+3],c[j+3],r3);
  }
}

void GF_transposed_matmul_add_N_A_v2(GF c[AIMER_N], const GF a[AIMER_N],const uint64_t bl[AIM2_NUM_BITS_FIELD],const uint64_t bh[AIM2_NUM_BITS_FIELD])
{
  const __m512i one = _mm512_set1_epi64(0x01);

  uint8_t stateOut[16];
  
  __m512i Al,Ah, Tmp, Tmp2;
  __m128i aux2;
  
  for(int j=0;j<16;j++){
    aux2 = _mm_loadu_si128((const __m128i *)a[j]);
    Al =  _mm512_broadcast_i32x2(aux2);
    Ah =  _mm512_broadcast_i32x2(_mm_bsrli_si128(aux2,8));

    for(int i=15,k=0;i>=0;i--,k++){
      Tmp  = _mm512_loadu_si512((const __m512i *)bl+i);
      Tmp2 = _mm512_loadu_si512((const __m512i *)bh+i);
        
      Tmp = _mm512_xor_epi64(_mm512_and_epi64(Al,Tmp),_mm512_and_epi64(Ah,Tmp2));
      Tmp = _mm512_popcnt_epi64(Tmp);
      Tmp = _mm512_and_epi64(Tmp,one);
          uint8_t o = _mm512_cmp_epi64_mask (Tmp, one, 0x00);
          stateOut[k] = o;
    }
    GF r;
    GF_from_bytes(r,stateOut);
    GF_add(c[j],c[j],r);
  }
}



void GF_transposed_matmul_add_N_A_v1(GF c[AIMER_N], const GF a[AIMER_N],const uint64_t bl[AIM2_NUM_BITS_FIELD],const uint64_t bh[AIM2_NUM_BITS_FIELD])
{
  const __m512i one = _mm512_set1_epi64(0x01);

  uint8_t stateOut[16];
  
  __m512i Al,Ah, Aux, Tmp, Tmp2;
  __m128i aux2;
  
  for(int j=0;j<4;j++){
    Aux = _mm512_loadu_si512((const __m512i *)a[j*4]);
    //-----------
    aux2 = _mm512_extracti32x4_epi32(Aux,0);
    Al =  _mm512_broadcast_i32x2(aux2);
    Ah =  _mm512_broadcast_i32x2(_mm_bsrli_si128(aux2,8));

    for(int i=15,k=0;i>=0;i--,k++){
      Tmp  = _mm512_loadu_si512((const __m512i *)bl+i);
      Tmp2 = _mm512_loadu_si512((const __m512i *)bh+i);
        
      Tmp = _mm512_xor_epi64(_mm512_and_epi64(Al,Tmp),_mm512_and_epi64(Ah,Tmp2));
      Tmp = _mm512_popcnt_epi64(Tmp);
      Tmp = _mm512_and_epi64(Tmp,one);
          uint8_t o = _mm512_cmp_epi64_mask (Tmp, one, 0x00);
          stateOut[k] = o;
    }
    GF r;
    GF_from_bytes(r,stateOut);
    GF_add(c[j*4],c[j*4],r);
    //-----------
    aux2 = _mm512_extracti32x4_epi32(Aux,1);
    Al =  _mm512_broadcast_i32x2(aux2);
    Ah =  _mm512_broadcast_i32x2(_mm_bsrli_si128(aux2,8));

    for(int i=15,k=0;i>=0;i--,k++){
      Tmp  = _mm512_loadu_si512((const __m512i *)bl+i);
      Tmp2 = _mm512_loadu_si512((const __m512i *)bh+i);
        
      Tmp = _mm512_xor_epi64(_mm512_and_epi64(Al,Tmp),_mm512_and_epi64(Ah,Tmp2));
      Tmp = _mm512_popcnt_epi64(Tmp);
      Tmp = _mm512_and_epi64(Tmp,one);
          uint8_t o = _mm512_cmp_epi64_mask (Tmp, one, 0x00);
          stateOut[k] = o;
    }
    GF_from_bytes(r,stateOut);
    GF_add(c[j*4+1],c[j*4+1],r);
    //------------------------------
    aux2 = _mm512_extracti32x4_epi32(Aux,2);
    Al =  _mm512_broadcast_i32x2(aux2);
    Ah =  _mm512_broadcast_i32x2(_mm_bsrli_si128(aux2,8));

    for(int i=15,k=0;i>=0;i--,k++){
      Tmp  = _mm512_loadu_si512((const __m512i *)bl+i);
      Tmp2 = _mm512_loadu_si512((const __m512i *)bh+i);
        
      Tmp = _mm512_xor_epi64(_mm512_and_epi64(Al,Tmp),_mm512_and_epi64(Ah,Tmp2));
      Tmp = _mm512_popcnt_epi64(Tmp);
      Tmp = _mm512_and_epi64(Tmp,one);
          uint8_t o = _mm512_cmp_epi64_mask (Tmp, one, 0x00);
          stateOut[k] = o;
    }
    GF_from_bytes(r,stateOut);
    GF_add(c[j*4+2],c[j*4+2],r);
    //-------------------------------
    aux2 = _mm512_extracti32x4_epi32(Aux,3);
    Al =  _mm512_broadcast_i32x2(aux2);
    Ah =  _mm512_broadcast_i32x2(_mm_bsrli_si128(aux2,8));

    for(int i=15,k=0;i>=0;i--,k++){
      Tmp  = _mm512_loadu_si512((const __m512i *)bl+i);
      Tmp2 = _mm512_loadu_si512((const __m512i *)bh+i);
        
      Tmp = _mm512_xor_epi64(_mm512_and_epi64(Al,Tmp),_mm512_and_epi64(Ah,Tmp2));
      Tmp = _mm512_popcnt_epi64(Tmp);
      Tmp = _mm512_and_epi64(Tmp,one);
          uint8_t o = _mm512_cmp_epi64_mask (Tmp, one, 0x00);
          stateOut[k] = o;
    }
    GF_from_bytes(r,stateOut);
    GF_add(c[j*4+3],c[j*4+3],r);
    //--------------------
  }

}

// static uint8_t bitrev(uint8_t b)
// {
//   uint64_t r = 0;
//   unsigned byte_len = 8;

//   while (byte_len--)
//   {
//     r = (r << 1) | (b & 1);
//     b >>= 1;
//   }
//   return r;
// }




// void GF_transposed_matmul_add_N_M_old(GF c[AIMER_N], const GF a[AIMER_N],const GF b[AIM2_NUM_BITS_FIELD])
// {
//   const __m256i zero = _mm256_setzero_si256();
//   __m256i m0, m1, m2, m3, t;
//   const __m512i one = _mm512_set1_epi64(0x01);

//   uint8_t stateOut[16];  
  

//   __m512i Al,Ah, Aux, Tmp, Tmp2;
//   __m128i aux2;
  
//   for(int j=0;j<4;j++){
//     Aux = _mm512_loadu_si512((const __m512i *)a[j*4]);
//     print_512(Aux);
//     //-----------
//     aux2 = _mm512_extracti32x4_epi32(Aux,0);
//     Al =  _mm512_broadcast_i32x2(aux2);
//     Ah =  _mm512_broadcast_i32x2(_mm_bsrli_si128(aux2,8));
//     print_512(Al);
//     print_512(Ah);
//     getchar();

//     for(int i=0;i<16;i++){
//       Tmp  = _mm512_loadu_si512((const __m512i *)b[i]);
//       Tmp2 = _mm512_loadu_si512((const __m512i *)b[64+i]);
        
//       Tmp = _mm512_xor_epi64(_mm512_and_epi64(Al,Tmp),_mm512_and_epi64(Ah,Tmp2));
//       Tmp = _mm512_popcnt_epi64(Tmp);
//       Tmp = _mm512_and_epi64(Tmp,one);
//           uint8_t o = _mm512_cmp_epi64_mask (Tmp, one, 0x00);
//           stateOut[i] = o;
//           printf("%.2x\n",o);
//           getchar();
//     }
//     GF_add(c[j*4],c[j*4],stateOut);
//     //-----------
//     aux2 = _mm512_extracti32x4_epi32(Aux,1);
//     Al =  _mm512_broadcast_i32x2(aux2);
//     Ah =  _mm512_broadcast_i32x2(_mm_bsrli_si128(aux2,8));

//     for(int i=0;i<16;i++){
//       Tmp  = _mm512_loadu_si512((const __m512i *)b[i]);
//       Tmp2 = _mm512_loadu_si512((const __m512i *)b[64+i]);
        
//       Tmp = _mm512_xor_epi64(_mm512_and_epi64(Al,Tmp),_mm512_and_epi64(Ah,Tmp2));
//       Tmp = _mm512_popcnt_epi64(Tmp);
//       Tmp = _mm512_and_epi64(Tmp,one);
//           uint8_t o = _mm512_cmp_epi64_mask (Tmp, one, 0x00);
//           stateOut[i] = o;
//     }
//     GF_add(c[j*4+1],c[j*4+1],stateOut);
//     //------------------------------
//     aux2 = _mm512_extracti32x4_epi32(Aux,2);
//     Al =  _mm512_broadcast_i32x2(aux2);
//     Ah =  _mm512_broadcast_i32x2(_mm_bsrli_si128(aux2,8));

//     for(int i=0;i<16;i++){
//       Tmp  = _mm512_loadu_si512((const __m512i *)b[i]);
//       Tmp2 = _mm512_loadu_si512((const __m512i *)b[64+i]);
        
//       Tmp = _mm512_xor_epi64(_mm512_and_epi64(Al,Tmp),_mm512_and_epi64(Ah,Tmp2));
//       Tmp = _mm512_popcnt_epi64(Tmp);
//       Tmp = _mm512_and_epi64(Tmp,one);
//           uint8_t o = _mm512_cmp_epi64_mask (Tmp, one, 0x00);
//           stateOut[i] = o;
//     }
//     GF_add(c[j*4+2],c[j*4+2],stateOut);
//     //-------------------------------
//     aux2 = _mm512_extracti32x4_epi32(Aux,3);
//     Al =  _mm512_broadcast_i32x2(aux2);
//     Ah =  _mm512_broadcast_i32x2(_mm_bsrli_si128(aux2,8));

//     for(int i=0;i<16;i++){
//       Tmp  = _mm512_loadu_si512((const __m512i *)b[i]);
//       Tmp2 = _mm512_loadu_si512((const __m512i *)b[64+i]);
        
//       Tmp = _mm512_xor_epi64(_mm512_and_epi64(Al,Tmp),_mm512_and_epi64(Ah,Tmp2));
//       Tmp = _mm512_popcnt_epi64(Tmp);
//       Tmp = _mm512_and_epi64(Tmp,one);
//           uint8_t o = _mm512_cmp_epi64_mask (Tmp, one, 0x00);
//           stateOut[i] = o;
//     }
//     GF_add(c[j*4+3],c[j*4+3],stateOut);
//     //--------------------
//   }


//   // for(int i=0;i<16;i++)  
//   //   printBytes(c[i],16);
    
//   // getchar();

// }

void GF_transposed_matmul_add_N(GF c[AIMER_N], const GF a[AIMER_N],
                                const GF b[AIM2_NUM_BITS_FIELD])
{
  const __m256i zero = _mm256_setzero_si256();
  __m256i m0, m1, m2, m3, t;

  // printBytes(c[0],16);
  // printBytes(a[1],16);
  // printBytes(a[2],16);
  // printBytes(a[3],16);
  // printBytes(b[0],16);
  // getchar();

  for (size_t party = 0; party < AIMER_N; party += 4)
  {
    const __m256i shift = _mm256_set_epi64x(0, 0, 1, 1);
    __m256i a0, a1, a2, a3, c0, c1, c2, c3;

    m0 = _mm256_loadu_si256((const __m256i *)a[party]);
    m1 = _mm256_loadu_si256((const __m256i *)a[party + 2]);
    // m2 = _mm256_setzero_si256();
    // m3 = _mm256_setzero_si256();
    m2 = _mm256_loadu_si256((const __m256i *)c[party]);
    m3 = _mm256_loadu_si256((const __m256i *)c[party + 2]);

    // printf("party=%d\n",party);
    // print_256(m0);
    // print_256(m2);
    // printBytes(a[0],16);
    // printBytes(c[0],16);
    // getchar();

    // ai = a[party + i].hi x4
    a0 = _mm256_permute4x64_epi64(m0, 0x55);
    a1 = _mm256_permute4x64_epi64(m0, 0xff);
    a2 = _mm256_permute4x64_epi64(m1, 0x55);
    a3 = _mm256_permute4x64_epi64(m1, 0xff);

    // print_256(a0);
    // print_256(a1);
    // print_256(a2);
    // print_256(a3);
    // getchar();

    a0 = _mm256_sllv_epi64(a0, shift);
    a1 = _mm256_sllv_epi64(a1, shift);
    a2 = _mm256_sllv_epi64(a2, shift);
    a3 = _mm256_sllv_epi64(a3, shift);

    // print_256(a0);
    // print_256(a1);
    // print_256(a2);
    // print_256(a3);
    // getchar();

    // ci = [c[party + i].lo | c[party + i].hi | 0 | 0]
    c0 = _mm256_permute2x128_si256(m2, m2, 0x80);
    c1 = _mm256_permute2x128_si256(m2, m2, 0x81);
    c2 = _mm256_permute2x128_si256(m3, m3, 0x80);
    c3 = _mm256_permute2x128_si256(m3, m3, 0x81);
    // print_256(c0);
    // getchar();

    for (size_t row = 128; row > 64; row -= 2)
    {
      t = _mm256_loadu_si256((const __m256i *)b[row - 2]);

      // mi[w] = (ai[w].msb) ^ 64
      m0 = _mm256_cmpgt_epi64(zero, a0);
      m1 = _mm256_cmpgt_epi64(zero, a1);
      m2 = _mm256_cmpgt_epi64(zero, a2);
      m3 = _mm256_cmpgt_epi64(zero, a3);

      // print_256(m0);
      // print_256(m1);
      

      a0 = _mm256_slli_epi64(a0, 2);
      a1 = _mm256_slli_epi64(a1, 2);
      a2 = _mm256_slli_epi64(a2, 2);
      a3 = _mm256_slli_epi64(a3, 2);

      // print_256(t);
      // print_256(m1);
      // getchar();

      m0 = _mm256_and_si256(m0, t);
      m1 = _mm256_and_si256(m1, t);
      m2 = _mm256_and_si256(m2, t);
      m3 = _mm256_and_si256(m3, t);

      // print_256(m0);
      

      c0 = _mm256_xor_si256(c0, m0);
      c1 = _mm256_xor_si256(c1, m1);
      c2 = _mm256_xor_si256(c2, m2);
      c3 = _mm256_xor_si256(c3, m3);
    }
    // getchar();

    m0 = _mm256_loadu_si256((const __m256i *)a[party]);
    m1 = _mm256_loadu_si256((const __m256i *)a[party + 2]);

    // ai = a[party + i].lo x4
    a0 = _mm256_permute4x64_epi64(m0, 0x00);
    a1 = _mm256_permute4x64_epi64(m0, 0xaa);
    a2 = _mm256_permute4x64_epi64(m1, 0x00);
    a3 = _mm256_permute4x64_epi64(m1, 0xaa);

    a0 = _mm256_sllv_epi64(a0, shift);
    a1 = _mm256_sllv_epi64(a1, shift);
    a2 = _mm256_sllv_epi64(a2, shift);
    a3 = _mm256_sllv_epi64(a3, shift);

    for (size_t row = 64; row > 0; row -= 2)
    {
      t = _mm256_loadu_si256((const __m256i *)b[row - 2]);

      // mi[w] = (ai[w].msb) ^ 64
      m0 = _mm256_cmpgt_epi64(zero, a0);
      m1 = _mm256_cmpgt_epi64(zero, a1);
      m2 = _mm256_cmpgt_epi64(zero, a2);
      m3 = _mm256_cmpgt_epi64(zero, a3);

      // print_256(m0);
      // print_256(m1);

      a0 = _mm256_slli_epi64(a0, 2);
      a1 = _mm256_slli_epi64(a1, 2);
      a2 = _mm256_slli_epi64(a2, 2);
      a3 = _mm256_slli_epi64(a3, 2);

      

      m0 = _mm256_and_si256(m0, t);
      m1 = _mm256_and_si256(m1, t);
      m2 = _mm256_and_si256(m2, t);
      m3 = _mm256_and_si256(m3, t);
      
      // print_256(m0);
      

      c0 = _mm256_xor_si256(c0, m0);
      c1 = _mm256_xor_si256(c1, m1);
      c2 = _mm256_xor_si256(c2, m2);
      c3 = _mm256_xor_si256(c3, m3);
    }
    // getchar();

    a0 = _mm256_permute2x128_si256(c0, c1, 0x20);
    a1 = _mm256_permute2x128_si256(c0, c1, 0x31);
    a2 = _mm256_permute2x128_si256(c2, c3, 0x20);
    a3 = _mm256_permute2x128_si256(c2, c3, 0x31);

    c0 = _mm256_xor_si256(a0, a1);
    c2 = _mm256_xor_si256(a2, a3);
    // print_256(c0);
    // getchar();

    _mm256_storeu_si256((__m256i *)c[party], c0);
    _mm256_storeu_si256((__m256i *)c[party + 2], c2);
  }
  // printBytes(c[0],16);
  // printBytes(c[1],16);
  // printBytes(c[2],16);
  // printBytes(c[3],16);
  // getchar();
}

void POLY_mul_add_N(GF lo[AIMER_N], GF hi[AIMER_N],
                    const GF a[AIMER_N], const GF b)
{
  __m128i x[2], y, z[2], zhi[2], t[4];
  y = _mm_loadu_si128((const __m128i *)b);

  for (size_t party = 0; party < AIMER_N; party += 2)
  {
    x[0] = _mm_loadu_si128((const __m128i *)a[party]);
    x[1] = _mm_loadu_si128((const __m128i *)a[party + 1]);

    // polynomial multiplication x2
    t[0] = _mm_clmulepi64_si128(x[0], y, 0x01);
    t[1] = _mm_clmulepi64_si128(x[1], y, 0x01);
    t[2] = _mm_clmulepi64_si128(x[0], y, 0x10);
    t[3] = _mm_clmulepi64_si128(x[1], y, 0x10);

    z[0] = _mm_clmulepi64_si128(x[0], y, 0x00);
    z[1] = _mm_clmulepi64_si128(x[1], y, 0x00);
    zhi[0] = _mm_clmulepi64_si128(x[0], y, 0x11);
    zhi[1] = _mm_clmulepi64_si128(x[1], y, 0x11);

    z[0] = _mm_xor_si128(*(__m128i *)lo[party], z[0]);
    z[1] = _mm_xor_si128(*(__m128i *)lo[party + 1], z[1]);
    zhi[0] = _mm_xor_si128(*(__m128i *)hi[party], zhi[0]);
    zhi[1] = _mm_xor_si128(*(__m128i *)hi[party + 1], zhi[1]);

    t[0] = _mm_xor_si128(t[0], t[2]);
    t[1] = _mm_xor_si128(t[1], t[3]);

    z[0] = _mm_xor_si128(z[0], _mm_slli_si128(t[0], 8));
    z[1] = _mm_xor_si128(z[1], _mm_slli_si128(t[1], 8));
    zhi[0] = _mm_xor_si128(zhi[0], _mm_srli_si128(t[0], 8));
    zhi[1] = _mm_xor_si128(zhi[1], _mm_srli_si128(t[1], 8));

    _mm_storeu_si128((__m128i *)lo[party], z[0]);
    _mm_storeu_si128((__m128i *)lo[party + 1], z[1]);
    _mm_storeu_si128((__m128i *)hi[party], zhi[0]);
    _mm_storeu_si128((__m128i *)hi[party + 1], zhi[1]);
  }
}

void POLY_red_N(GF lo[AIMER_N], const GF hi[AIMER_N])
{
  __m128i z[4], zhi[4], t[4];
  __m128i irr = _mm_set_epi64x(0x0, 0x87);

  for (size_t party = 0; party < AIMER_N; party += 4)
  {
    z[0] = _mm_loadu_si128((const __m128i *)lo[party]);
    z[1] = _mm_loadu_si128((const __m128i *)lo[party + 1]);
    z[2] = _mm_loadu_si128((const __m128i *)lo[party + 2]);
    z[3] = _mm_loadu_si128((const __m128i *)lo[party + 3]);

    zhi[0] = _mm_loadu_si128((const __m128i *)hi[party]);
    zhi[1] = _mm_loadu_si128((const __m128i *)hi[party + 1]);
    zhi[2] = _mm_loadu_si128((const __m128i *)hi[party + 2]);
    zhi[3] = _mm_loadu_si128((const __m128i *)hi[party + 3]);

    // modular reduction x4
    t[0] = _mm_clmulepi64_si128(zhi[0], irr, 0x01);
    t[1] = _mm_clmulepi64_si128(zhi[1], irr, 0x01);
    t[2] = _mm_clmulepi64_si128(zhi[2], irr, 0x01);
    t[3] = _mm_clmulepi64_si128(zhi[3], irr, 0x01);

    z[0] = _mm_xor_si128(z[0], _mm_slli_si128(t[0], 8));
    z[1] = _mm_xor_si128(z[1], _mm_slli_si128(t[1], 8));
    zhi[0] = _mm_xor_si128(zhi[0], _mm_srli_si128(t[0], 8));
    zhi[1] = _mm_xor_si128(zhi[1], _mm_srli_si128(t[1], 8));

    z[2] = _mm_xor_si128(z[2], _mm_slli_si128(t[2], 8));
    z[3] = _mm_xor_si128(z[3], _mm_slli_si128(t[3], 8));
    zhi[2] = _mm_xor_si128(zhi[2], _mm_srli_si128(t[2], 8));
    zhi[3] = _mm_xor_si128(zhi[3], _mm_srli_si128(t[3], 8));

    t[0] = _mm_clmulepi64_si128(zhi[0], irr, 0x00);
    t[1] = _mm_clmulepi64_si128(zhi[1], irr, 0x00);
    t[2] = _mm_clmulepi64_si128(zhi[2], irr, 0x00);
    t[3] = _mm_clmulepi64_si128(zhi[3], irr, 0x00);

    z[0] = _mm_xor_si128(z[0], t[0]);
    z[1] = _mm_xor_si128(z[1], t[1]);
    z[2] = _mm_xor_si128(z[2], t[2]);
    z[3] = _mm_xor_si128(z[3], t[3]);

    _mm_storeu_si128((__m128i *)lo[party], z[0]);
    _mm_storeu_si128((__m128i *)lo[party + 1], z[1]);
    _mm_storeu_si128((__m128i *)lo[party + 2], z[2]);
    _mm_storeu_si128((__m128i *)lo[party + 3], z[3]);
  }
}
