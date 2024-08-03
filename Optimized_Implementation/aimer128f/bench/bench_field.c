#include "bench.h"
#include "../field.h"
#include "../aim2.h"
#include "../common/rng.h"

static void GF_rand(GF a)
{
    uint8_t random[AIM2_NUM_BYTES_FIELD];
    randombytes(random, AIM2_NUM_BYTES_FIELD);
    GF_from_bytes(a, random);
}

static void GF_rand_N(GF a[AIMER_N])
{
    for (int i = 0; i < AIMER_N; i++)
    {
        GF_rand(a[i]);
    }
}

void bench_field(void)
{
    GF a, b, c;
    GF A[AIMER_N], C[AIMER_N];
    GF LO[AIMER_N], HI[AIMER_N];
    GF matrix_A[AIMER_L][AIM2_NUM_BITS_FIELD];
    GF vector_b = {0};
    uint8_t iv[AIM2_IV_SIZE];

    GF_rand(a);
    GF_rand(b);
    GF_rand(c);
    GF_rand_N(A);
    GF_rand_N(C);
    generate_matrix_LU(matrix_A, vector_b, iv);

    int BENCH = 3000;
    BENCH_FUNCTION(GF_add, c, a, b);
    BENCH_FUNCTION(GF_mul, c, a, b);
    // New functions
    BENCH_FUNCTION(GF_mul_red, c, a, b);
    BENCH_FUNCTION(GF_mul_karat, c, a, b);
    BENCH_FUNCTION(GF_mul_karat_red, c, a, b);

    BENCH_FUNCTION(GF_sqr, c, a);
    BENCH_FUNCTION(GF_mul_N, C, A, b);
   
    // New function
    BENCH_FUNCTION(GF_mul_N_karat_red, C, A, b);

    BENCH_FUNCTION(GF_sqr_N, C, A);
    BENCH_FUNCTION(POLY_mul_add_N, LO, HI, A, b);
    BENCH_FUNCTION(POLY_red_N, LO, HI);
    BENCH = 500;
    BENCH_FUNCTION(GF_transposed_matmul, c, a, (const GF(*))matrix_A);
    BENCH_FUNCTION(GF_transposed_matmul_add_N, C, A, (const GF(*))matrix_A);
}
