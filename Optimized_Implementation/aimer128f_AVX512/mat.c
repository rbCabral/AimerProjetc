// SPDX-License-Identifier: MIT

#include <stdio.h>
#include <string.h>
#include "field.h"
#include "aim2.h"

uint64_t bitrev(uint64_t b)
{
  uint64_t r = 0;
  unsigned byte_len = 64;

  while (byte_len--)
  {
    r = (r << 1) | (b & 1);
    b >>= 1;
  }
  return r;
}

void hex(char *s, GF x)
{
  printf("%s: %016lx,%016lx\n", s, x[1], x[0]);
}

void bin(GF x)
{
  printf("\"%064lb%064lb\",\n", bitrev(x[0]), bitrev(x[1]));
}

void mat(char *s, GF matrix[AIM2_NUM_BITS_FIELD])
{
  printf("%s = [\n", s);
  for (int i = 0; i < AIM2_NUM_BITS_FIELD; i++)
  {
    bin(matrix[i]);
  }
  printf("]\n");
}

int main()
{
  GF matrix_L[AIM2_NUM_INPUT_SBOX][AIM2_NUM_BITS_FIELD];
  GF matrix_U[AIM2_NUM_INPUT_SBOX][AIM2_NUM_BITS_FIELD];
  GF vector_b = {
      0,
  };
  const uint8_t iv[AIM2_IV_SIZE] = {0};
  GF matrix_A[AIM2_NUM_INPUT_SBOX][AIM2_NUM_BITS_FIELD];

  generate_matrices_L_and_U2(matrix_L, matrix_U, vector_b, iv);

  mat("LT", matrix_L[0]);
  mat("UT", matrix_U[0]);

  for (size_t i = 0; i < AIM2_NUM_BITS_FIELD; i++)
  {
    GF_transposed_matmul(matrix_A[0][i], matrix_U[0][i], (const GF *)matrix_L[0]);
  }

  mat("AT", matrix_A[0]);

  return 0;
}
