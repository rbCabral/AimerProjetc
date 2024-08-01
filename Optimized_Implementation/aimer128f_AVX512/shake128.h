#include <immintrin.h>

#define ALIGN2 __attribute__ ((aligned (32)))

void keccakf(__m512i *s);
// void print_512(__m512i x);
