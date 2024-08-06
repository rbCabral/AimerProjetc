// SPDX-License-Identifier: MIT

#include "field.h"
#include <stddef.h>
#include <stdint.h>
#include <immintrin.h>

#include <stdio.h>

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

static void printBytes(const uint8_t* in, int len){
    int i;
    for(i=0;i<len;i++){
        printf("%.2X", in[i]);
    }
    printf("\n");
}

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



void GF_transposed_matmul2(GF c, const GF a, const GF b[AIM2_NUM_BITS_FIELD])
{
  unsigned int i, j;
  const uint32_t *a_ptr = (uint32_t *)a;
  const __m256i shift = _mm256_set_epi32(0, 2, 4, 6, 1, 3, 5, 7);
  __m256i temp[2] = {_mm256_setzero_si256(), _mm256_setzero_si256()};
  __m256i matrix_data[4];


  for (i = 0; i < 4; i++)
  {
    __m256i index = _mm256_set1_epi32(a_ptr[i]);
    print_256(index);
    for (j = 32 * (i + 1); j > 32 * i; j -= 8)
    {
      __m256i mask = _mm256_sllv_epi32(index, shift);
      printf("id = %d\n",j - 2);
      matrix_data[0] = _mm256_loadu_si256((const __m256i *)&b[j - 2]);
      matrix_data[1] = _mm256_loadu_si256((const __m256i *)&b[j - 4]);
      matrix_data[2] = _mm256_loadu_si256((const __m256i *)&b[j - 6]);
      matrix_data[3] = _mm256_loadu_si256((const __m256i *)&b[j - 8]);

      print_256(matrix_data[0]);
      print_256(_mm256_srai_epi32(_mm256_shuffle_epi32(mask, 0xff), 31));
      print_256(matrix_data[1]);
      print_256(_mm256_srai_epi32(_mm256_shuffle_epi32(mask, 0xaa), 31));
      print_256(matrix_data[2]);
      print_256(_mm256_srai_epi32(_mm256_shuffle_epi32(mask, 0x55), 31));
      print_256(matrix_data[3]);
      print_256(_mm256_srai_epi32(_mm256_shuffle_epi32(mask, 0x00), 31));
      getchar();

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
  printf("fim\n");
  getchar();
  _mm_storeu_si128((__m128i *)c,
                  _mm_xor_si128(_mm256_extracti128_si256(temp[0], 0),
                  _mm256_extracti128_si256(temp[0], 1)));
  printBytes(c,16);
  getchar();                  
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

void GF_transposed_matmul_2(GF c, const GF a, const GF b[AIM2_NUM_BITS_FIELD])
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

// void GF_transposed_matmul_add_N_A(GF c[AIMER_N], const GF a[AIMER_N],const uint64_t b[AIM2_NUM_BITS_FIELD])
// {
//   const __m256i zero = _mm256_setzero_si256();
//   __m256i m0, m1, m2, m3, t;
//   uint8_t sLow[64];
//   uint8_t sHigh[64];

//   __m512i A1,A2, Ac1, Ac2;

//   // __m512i A = _mm_set_epi8(0xff,0xfc,0xf3,0xf0,0xcf,0xcc,0xc3,0xc0,0x3f,0x3c,0x33,0x30,0x0f,0x0c,0x03,0x00);
//   __m512i Mask = _mm512_set_epi64(0xfffcf3f0cfccc3c0,0x3f3c33300f0c0300,0xfffcf3f0cfccc3c0,0x3f3c33300f0c0300,0xfffcf3f0cfccc3c0,0x3f3c33300f0c0300,0xfffcf3f0cfccc3c0,0x3f3c33300f0c0300);
// 	// __m512i B = _mm_set_epi8(0xaa,0x10,0x12,0x9f,0x19,0xf5,0x9d,0x85,0x33,0x11,0xc1,0xb2,0xf3,0x14,0xe5,0x6a);
// 	__m512i M = _mm512_set_epi64(0x0707060605050404,0x0303020201010000,0x0707060605050404,0x0303020201010000,0x0707060605050404,0x0303020201010000,0x0707060605050404,0x0303020201010000);
// 	__m512i Tmp,Ok1,Ok2;

//   __m512i ZeroA =  _mm512_set1_epi16(0x000F);
// 	__m512i Zero2 =  _mm512_set1_epi16(0x0F00);

//   A1 = _mm512_loadu_si512((const __m512i *)a[0]);
//   A2 = _mm512_bsrli_epi128(A1,8);

//   print_512(A1);


// 	A1 = _mm512_shuffle_epi8(A1,M);
// 	Ok1 = _mm512_and_si512(A1,ZeroA);
// 	Ok1 = _mm512_shuffle_epi8(Mask,Ok1);

// 	Tmp = _mm512_srli_epi16(A1,4);
// 	Ok2 = _mm512_and_si512(Tmp,Zero2);
// 	Ok2 = _mm512_shuffle_epi8(Mask,Ok2);
// 	A1  = _mm512_xor_si512(Ok1,Ok2);

//   A2 = _mm512_shuffle_epi8(A2,M);
// 	Ok1 = _mm512_and_si512(A2,ZeroA);
// 	Ok1 = _mm512_shuffle_epi8(Mask,Ok1);

// 	Tmp = _mm512_srli_epi16(A2,4);
// 	Ok2 = _mm512_and_si512(Tmp,Zero2);
// 	Ok2 = _mm512_shuffle_epi8(Mask,Ok2);
// 	A2 = _mm512_xor_si512(Ok1,Ok2);
	
//   // print_512(A1);
//   // print_512(A2);
//   // getchar();
//   _mm512_storeu_si512((__m512i *)sLow, A1);
//   _mm512_storeu_si512((__m512i *)sHigh, A2);
//   Ac1 = _mm512_setzero_si512();
//   Ac2 = _mm512_setzero_si512();

  
//   // print_512(A1);
//   // print_512(A2);

//   // printf("%.2x\n",sLow[0]);
//   // printf("%.2x\n",sLow[1]);
//   // printf("%.2x\n",sLow[2]);
//   // printf("%.2x\n",sLow[3]);
//   // printf("%.2x\n",sHigh[0]);
  
//   printf("aqui\n");

//   for(int i=0, j=0;i<64;i++,j++){
//     print_512(_mm512_maskz_loadu_epi64(sLow[i] ,(const __m512i *)b+j));
//     // print_512(_mm512_maskz_loadu_epi64(sHigh[i] ,(const __m512i *)b+j+1));
//     Ac1 = _mm512_xor_si512(Ac1,_mm512_maskz_loadu_epi64(sLow[i] ,(const __m512i *)b+j));
//     Ac2 = _mm512_xor_si512(Ac2,_mm512_maskz_loadu_epi64(sHigh[i],(const __m512i *)b+j+1));
//   }
//   getchar();

//   A1 = _mm512_loadu_si512((const __m512i *)c);
//   Ok1 = _mm512_xor_si512(Ac1,Ac2);
//   // print_512(Ok1);
//   Ok1 = _mm512_xor_si512(Ok1,A1);
//   // print_512(Ok1);
//   // getchar();

//   _mm512_storeu_si512((__m512i *)c[0], Ok1);
//   printBytes(c[0],16);
//   printBytes(c[1],16);
//   printBytes(c[2],16);
//   printBytes(c[3],16);
//   getchar();

  
  

void GF_transposed_matmul_add_N_A(GF c[AIMER_N], const GF a[AIMER_N],const uint64_t bl[AIM2_NUM_BITS_FIELD],const uint64_t bh[AIM2_NUM_BITS_FIELD])
{
  const __m256i zero = _mm256_setzero_si256();
  __m256i m0, m1, m2, m3, t;
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
    GF_add(c[j*4],c[j*4],stateOut);
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
    GF_add(c[j*4+1],c[j*4+1],stateOut);
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
    GF_add(c[j*4+2],c[j*4+2],stateOut);
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
    GF_add(c[j*4+3],c[j*4+3],stateOut);
    //--------------------
  }


  // for(int i=0;i<16;i++)  
  //   printBytes(c[i],16);
    
  // getchar();

}


void GF_transposed_matmul_add_N_M(GF c[AIMER_N], const GF a[AIMER_N],const GF b[AIM2_NUM_BITS_FIELD])
{
  const __m256i zero = _mm256_setzero_si256();
  __m256i m0, m1, m2, m3, t;
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

    for(int i=0;i<16;i++){
      Tmp  = _mm512_loadu_si512((const __m512i *)b[i]);
      Tmp2 = _mm512_loadu_si512((const __m512i *)b[64+i]);
        
      Tmp = _mm512_xor_epi64(_mm512_and_epi64(Al,Tmp),_mm512_and_epi64(Ah,Tmp2));
      Tmp = _mm512_popcnt_epi64(Tmp);
      Tmp = _mm512_and_epi64(Tmp,one);
          uint8_t o = _mm512_cmp_epi64_mask (Tmp, one, 0x00);
          stateOut[i] = o;
    }
    GF_add(c[j*4],c[j*4],stateOut);
    //-----------
    aux2 = _mm512_extracti32x4_epi32(Aux,1);
    Al =  _mm512_broadcast_i32x2(aux2);
    Ah =  _mm512_broadcast_i32x2(_mm_bsrli_si128(aux2,8));

    for(int i=0;i<16;i++){
      Tmp  = _mm512_loadu_si512((const __m512i *)b[i]);
      Tmp2 = _mm512_loadu_si512((const __m512i *)b[64+i]);
        
      Tmp = _mm512_xor_epi64(_mm512_and_epi64(Al,Tmp),_mm512_and_epi64(Ah,Tmp2));
      Tmp = _mm512_popcnt_epi64(Tmp);
      Tmp = _mm512_and_epi64(Tmp,one);
          uint8_t o = _mm512_cmp_epi64_mask (Tmp, one, 0x00);
          stateOut[i] = o;
    }
    GF_add(c[j*4+1],c[j*4+1],stateOut);
    //------------------------------
    aux2 = _mm512_extracti32x4_epi32(Aux,2);
    Al =  _mm512_broadcast_i32x2(aux2);
    Ah =  _mm512_broadcast_i32x2(_mm_bsrli_si128(aux2,8));

    for(int i=0;i<16;i++){
      Tmp  = _mm512_loadu_si512((const __m512i *)b[i]);
      Tmp2 = _mm512_loadu_si512((const __m512i *)b[64+i]);
        
      Tmp = _mm512_xor_epi64(_mm512_and_epi64(Al,Tmp),_mm512_and_epi64(Ah,Tmp2));
      Tmp = _mm512_popcnt_epi64(Tmp);
      Tmp = _mm512_and_epi64(Tmp,one);
          uint8_t o = _mm512_cmp_epi64_mask (Tmp, one, 0x00);
          stateOut[i] = o;
    }
    GF_add(c[j*4+2],c[j*4+2],stateOut);
    //-------------------------------
    aux2 = _mm512_extracti32x4_epi32(Aux,3);
    Al =  _mm512_broadcast_i32x2(aux2);
    Ah =  _mm512_broadcast_i32x2(_mm_bsrli_si128(aux2,8));

    for(int i=0;i<16;i++){
      Tmp  = _mm512_loadu_si512((const __m512i *)b[i]);
      Tmp2 = _mm512_loadu_si512((const __m512i *)b[64+i]);
        
      Tmp = _mm512_xor_epi64(_mm512_and_epi64(Al,Tmp),_mm512_and_epi64(Ah,Tmp2));
      Tmp = _mm512_popcnt_epi64(Tmp);
      Tmp = _mm512_and_epi64(Tmp,one);
          uint8_t o = _mm512_cmp_epi64_mask (Tmp, one, 0x00);
          stateOut[i] = o;
    }
    GF_add(c[j*4+3],c[j*4+3],stateOut);
    //--------------------
  }


  // for(int i=0;i<16;i++)  
  //   printBytes(c[i],16);
    
  // getchar();

}

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
