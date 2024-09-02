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

void vet(char *s, GF matrix)
{
  printf("%s = [\n", s);
  bin(matrix);
  printf("]\n");
}