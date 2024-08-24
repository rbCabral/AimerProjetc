// SPDX-License-Identifier: MIT

#include "api.h"
#include "aim2.h"
#include "field.h"
#include "hash.h"
#include "params.h"
#include "sign.h"
#include "tree.h"
#include "common/crypto_declassify.h"
#include "common/rng.h"
#include <stddef.h>
#include <stdint.h>
#include <stdlib.h>
#include <string.h>

#include "bench.h"
#include "shake128.h"


#define BENCH 1500
// #define benchPhases

static void printBytes(const uint8_t* in, int len){
    int i;
    for(i=0;i<len;i++){
        printf("%.2X", in[i]);
    }
    printf("\n");
}

// void print_512(__m512i x){
//     // uint8_t c[64];
//     uint64_t c[8];
//     _mm512_storeu_si512((__m512i*)c+0,x);
//     printf("%8.16lx %8.16lx %8.16lx %8.16lx %8.16lx %8.16lx %8.16lx %8.16lx\n", c[7],c[6],c[5],c[4],c[3],c[2],c[1],c[0]);
//     // printBytes(c,64);
// }

void commit_and_expand_tape_x8(tape_t *tapes, uint8_t *commits,uint8_t* salt,const uint8_t *seeds, size_t rep, size_t party){
  __m512i s[25];
  // uint8_t bufs[4][AIMER_SEED_SIZE + 2];
  uint8_t bufs[4][64];
  uint8_t init[17];
  // uint8_t state[1600];
  // const uint8_t *in_ptrs[4] = {bufs[0], bufs[1], bufs[2], bufs[3]};
  const __m512i m0 = _mm512_set_epi64(0x00,0xFFFFFFFF,0xFFFFFFFFFFFFFFFF,0xFFFFFFFFFFFFFFFF,0x00,0xFFFFFFFF,0xFFFFFFFFFFFFFFFF,0xFFFFFFFFFFFFFFFF);
  
  init[0] = 0x05;
  memcpy(init+1, salt,16);
  s[0] = _mm512_set1_epi64(((uint64_t*)init)[0]);
  s[1] = _mm512_set1_epi64(((uint64_t*)init)[1]);
  
  for(int i=0;i<4;i++){
    bufs[i][0]  = init[16];//last byte of the salt  
    bufs[i][32] = init[16];  
    bufs[i][1]  = (uint8_t)(rep);//repetition
    bufs[i][33] = (uint8_t)(rep);
    bufs[i][2]  = (uint8_t)(party + i);//party
    bufs[i][34] = (uint8_t)(party + i + 4);
    memcpy(&bufs[i][3], seeds + (i)   * AIMER_SEED_SIZE, AIMER_SEED_SIZE);//seed
    memcpy(&bufs[i][35], seeds + (i+4) * AIMER_SEED_SIZE, AIMER_SEED_SIZE);
    bufs[i][19] = 0x1F;//shake128
    bufs[i][51] = 0x1F;  

  }

  s[2] = _mm512_and_si512(_mm512_loadu_si512((__m512i*)bufs[0]),m0);
  s[3] = _mm512_and_si512(_mm512_loadu_si512((__m512i*)bufs[1]),m0);
  s[4] = _mm512_and_si512(_mm512_loadu_si512((__m512i*)bufs[2]),m0);
  s[5] = _mm512_and_si512(_mm512_loadu_si512((__m512i*)bufs[3]),m0);

  s[6] = _mm512_unpackhi_epi64(s[2],s[3]);
  s[7] = _mm512_unpackhi_epi64(s[4],s[5]);
  s[2] = _mm512_unpacklo_epi64(s[2],s[3]);
  s[3] = _mm512_unpacklo_epi64(s[4],s[5]);
          
  s[4] = _mm512_permutex2var_epi64(s[2],_mm512_set_epi64(15,14,7,6,11,10,3,2),s[3]);
  s[5] = _mm512_permutex2var_epi64(s[6],_mm512_set_epi64(15,14,7,6,11,10,3,2),s[7]);
  s[2] = _mm512_permutex2var_epi64(s[2],_mm512_set_epi64(13,12,5,4,9,8,1,0),s[3]);
  s[3] = _mm512_permutex2var_epi64(s[6],_mm512_set_epi64(13,12,5,4,9,8,1,0),s[7]);

  for(int i=6;i<25;i++)
    s[i] = _mm512_setzero_si512();
  s[20] = _mm512_set1_epi64(0x8000000000000000);
  
  
  // for(int i=0;i<25;i++){
  //   _mm512_storeu_si512((__m512i*)state+i,s[i]);
  // }
  // printBytes(state,1600);
  keccakf(s);

  s[21] = _mm512_unpacklo_epi64(s[0],s[1]);
  s[22] = _mm512_unpackhi_epi64(s[0],s[1]);
	s[23] = _mm512_unpacklo_epi64(s[2],s[3]);
  s[24] = _mm512_unpackhi_epi64(s[2],s[3]);
  
  const __m512i l0 = _mm512_set_epi64(11,10,3,2,9,8,1,0);
  const __m512i l1 = _mm512_set_epi64(15,14,7,6,13,12,5,4);

	s[19] = _mm512_permutex2var_epi64(s[21],l0,s[23]);
	s[20] = _mm512_permutex2var_epi64(s[22],l0,s[24]);
	s[21] = _mm512_permutex2var_epi64(s[21],l1,s[23]);
	s[22] = _mm512_permutex2var_epi64(s[22],l1,s[24]);

  
  _mm256_storeu_si256((__m256i*)commits + 0,_mm512_extracti64x4_epi64(s[19],0));
  _mm256_storeu_si256((__m256i*)commits + 1,_mm512_extracti64x4_epi64(s[20],0));
  _mm256_storeu_si256((__m256i*)commits + 2,_mm512_extracti64x4_epi64(s[19],1));
  _mm256_storeu_si256((__m256i*)commits + 3,_mm512_extracti64x4_epi64(s[20],1));
  _mm256_storeu_si256((__m256i*)commits + 4,_mm512_extracti64x4_epi64(s[21],0));
  _mm256_storeu_si256((__m256i*)commits + 5,_mm512_extracti64x4_epi64(s[22],0));
  _mm256_storeu_si256((__m256i*)commits + 6,_mm512_extracti64x4_epi64(s[21],1));
  _mm256_storeu_si256((__m256i*)commits + 7,_mm512_extracti64x4_epi64(s[22],1));

  s[0] = _mm512_unpacklo_epi64(s[4],s[5]);
  s[1] = _mm512_unpackhi_epi64(s[4],s[5]);
	s[2] = _mm512_unpacklo_epi64(s[6],s[7]);
  s[3] = _mm512_unpackhi_epi64(s[6],s[7]);

  s[4] = _mm512_unpacklo_epi64(s[8],s[9]);
  s[5] = _mm512_unpackhi_epi64(s[8],s[9]);
	s[6] = _mm512_unpacklo_epi64(s[10],s[11]);
  s[7] = _mm512_unpackhi_epi64(s[10],s[11]);

  s[8] = _mm512_unpacklo_epi64(s[12],s[13]);
  s[9] = _mm512_unpackhi_epi64(s[12],s[13]);

  

  s[10] = _mm512_permutex2var_epi64(s[0],l0,s[2]);
	s[11] = _mm512_permutex2var_epi64(s[1],l0,s[3]);
	s[12] = _mm512_permutex2var_epi64(s[0],l1,s[2]);
	s[13] = _mm512_permutex2var_epi64(s[1],l1,s[3]);

  s[0] = _mm512_permutex2var_epi64(s[4],l0,s[6]);
	s[1] = _mm512_permutex2var_epi64(s[5],l0,s[7]);
	s[2] = _mm512_permutex2var_epi64(s[4],l1,s[6]);
	s[3] = _mm512_permutex2var_epi64(s[5],l1,s[7]);

  s[24] = _mm512_inserti64x4(s[10],_mm512_castsi512_si256(s[0]),1);
  s[23] = _mm512_inserti64x4(s[11],_mm512_castsi512_si256(s[1]),1);
  s[22] = _mm512_inserti64x4(s[0],_mm512_extracti64x4_epi64(s[10],1),0);
  s[21] = _mm512_inserti64x4(s[1],_mm512_extracti64x4_epi64(s[11],1),0);

  s[20] = _mm512_inserti64x4(s[12],_mm512_castsi512_si256(s[2]),1);
  s[19] = _mm512_inserti64x4(s[13],_mm512_castsi512_si256(s[3]),1);
  s[18] = _mm512_inserti64x4(s[2],_mm512_extracti64x4_epi64(s[12],1),0);
  s[17] = _mm512_inserti64x4(s[3],_mm512_extracti64x4_epi64(s[13],1),0);

  _mm512_storeu_si512((__m512i*)(uint8_t *)(tapes+0),s[24]);
  _mm512_storeu_si512((__m512i*)(uint8_t *)(tapes+1),s[23]);
  _mm512_storeu_si512((__m512i*)(uint8_t *)(tapes+2),s[22]);
  _mm512_storeu_si512((__m512i*)(uint8_t *)(tapes+3),s[21]);
  _mm512_storeu_si512((__m512i*)(uint8_t *)(tapes+4),s[20]);
  _mm512_storeu_si512((__m512i*)(uint8_t *)(tapes+5),s[19]);
  _mm512_storeu_si512((__m512i*)(uint8_t *)(tapes+6),s[18]);
  _mm512_storeu_si512((__m512i*)(uint8_t *)(tapes+7),s[17]);

  _mm_storeu_si128((__m128i*)(uint8_t *)(tapes+0)+4,_mm512_extracti64x2_epi64(s[8],0));
  _mm_storeu_si128((__m128i*)(uint8_t *)(tapes+1)+4,_mm512_extracti64x2_epi64(s[9],0));
  _mm_storeu_si128((__m128i*)(uint8_t *)(tapes+2)+4,_mm512_extracti64x2_epi64(s[8],1));
  _mm_storeu_si128((__m128i*)(uint8_t *)(tapes+3)+4,_mm512_extracti64x2_epi64(s[9],1));
  _mm_storeu_si128((__m128i*)(uint8_t *)(tapes+4)+4,_mm512_extracti64x2_epi64(s[8],2));
  _mm_storeu_si128((__m128i*)(uint8_t *)(tapes+5)+4,_mm512_extracti64x2_epi64(s[9],2));
  _mm_storeu_si128((__m128i*)(uint8_t *)(tapes+6)+4,_mm512_extracti64x2_epi64(s[8],3));
  _mm_storeu_si128((__m128i*)(uint8_t *)(tapes+7)+4,_mm512_extracti64x2_epi64(s[9],3));


  // printBytes((uint8_t *)(tapes),sizeof(tape_t));
  // printBytes((uint8_t *)(tapes+1),sizeof(tape_t));
  // printBytes((uint8_t *)(tapes+2),sizeof(tape_t));
  // printBytes((uint8_t *)(tapes+3),sizeof(tape_t));
  // printBytes((uint8_t *)(tapes+4),sizeof(tape_t));
  // printBytes((uint8_t *)(tapes+5),sizeof(tape_t));
  // printBytes((uint8_t *)(tapes+6),sizeof(tape_t));
  // printBytes((uint8_t *)(tapes+7),sizeof(tape_t));
  // getchar();

  // _mm512_storeu_si512((__m512i*)state,s[8]);
  // printBytes(state,64);
  // _mm512_storeu_si512((__m512i*)state,s[9]);
  // printBytes(state,64);
  // _mm512_storeu_si512((__m512i*)state,s[22]);
  // printBytes(state,64);
  // _mm512_storeu_si512((__m512i*)state,s[21]);
  // printBytes(state,64);
  
  // _mm512_storeu_si512((__m512i*)state,s[20]);
  // printBytes(state,64);
  // _mm512_storeu_si512((__m512i*)state,s[19]);
  // printBytes(state,64);
  // _mm512_storeu_si512((__m512i*)state,s[18]);
  // printBytes(state,64);
  // _mm512_storeu_si512((__m512i*)state,s[17]);
  // printBytes(state,64);
  

  


  // for(int i=10;i<14;i++){
  //   _mm512_storeu_si512((__m512i*)state,s[i]);
  //   printBytes(state,64);
  // }

  // for(int i=0;i<4;i++){
  //   _mm512_storeu_si512((__m512i*)state,s[i]);
  //   printBytes(state,64);
  // }
  
  // getchar();
	

  
  
  
  
  // getchar();

}


void commit_and_expand_tape_x4(tape_t *tapes, uint8_t *commits,
                               const hash_instance_x4 *ctx_precom,
                               const uint8_t *seeds, size_t rep, size_t party)
{
  hash_instance_x4 ctx;
  uint8_t bufs[4][AIMER_SEED_SIZE + 2];
  const uint8_t *in_ptrs[4] = {bufs[0], bufs[1], bufs[2], bufs[3]};
  uint8_t *out_ptrs[4];

  memcpy(&ctx, ctx_precom, sizeof(hash_instance_x4));

  bufs[0][0] = (uint8_t)(rep);
  bufs[1][0] = (uint8_t)(rep);
  bufs[2][0] = (uint8_t)(rep);
  bufs[3][0] = (uint8_t)(rep);

  bufs[0][1] = (uint8_t)(party + 0);
  bufs[1][1] = (uint8_t)(party + 1);
  bufs[2][1] = (uint8_t)(party + 2);
  bufs[3][1] = (uint8_t)(party + 3);

  memcpy(&bufs[0][2], seeds + 0 * AIMER_SEED_SIZE, AIMER_SEED_SIZE);
  memcpy(&bufs[1][2], seeds + 1 * AIMER_SEED_SIZE, AIMER_SEED_SIZE);
  memcpy(&bufs[2][2], seeds + 2 * AIMER_SEED_SIZE, AIMER_SEED_SIZE);
  memcpy(&bufs[3][2], seeds + 3 * AIMER_SEED_SIZE, AIMER_SEED_SIZE);

  // printBytes(bufs[0],1);
  // printBytes(bufs[1],1);
  // printBytes(bufs[2],1);
  // printBytes(bufs[3],1);
  // printBytes(&bufs[0][1],1);
  // printBytes(&bufs[1][1],1);
  // printBytes(&bufs[2][1],1);
  // printBytes(&bufs[3][1],1);
  // printBytes(&bufs[0][2],AIMER_SEED_SIZE);
  // printBytes(&bufs[1][2],AIMER_SEED_SIZE);
  // printBytes(&bufs[2][2],AIMER_SEED_SIZE);
  // printBytes(&bufs[3][2],AIMER_SEED_SIZE);
  // getchar();
  hash_update_x4(&ctx, in_ptrs, AIMER_SEED_SIZE + 2);
  hash_final_x4(&ctx);

  out_ptrs[0] = commits + 0 * AIMER_COMMIT_SIZE;
  out_ptrs[1] = commits + 1 * AIMER_COMMIT_SIZE;
  out_ptrs[2] = commits + 2 * AIMER_COMMIT_SIZE;
  out_ptrs[3] = commits + 3 * AIMER_COMMIT_SIZE;

  hash_squeeze_x4(&ctx, out_ptrs, AIMER_COMMIT_SIZE);
  
  // printBytes(out_ptrs[0],AIMER_COMMIT_SIZE);
  // printBytes(out_ptrs[1],AIMER_COMMIT_SIZE);
  // printBytes(out_ptrs[2],AIMER_COMMIT_SIZE);
  // printBytes(out_ptrs[3],AIMER_COMMIT_SIZE);

  out_ptrs[0] = (uint8_t *)(tapes);
  out_ptrs[1] = (uint8_t *)(tapes + 1);
  out_ptrs[2] = (uint8_t *)(tapes + 2);
  out_ptrs[3] = (uint8_t *)(tapes + 3);

  hash_squeeze_x4(&ctx, out_ptrs, sizeof(tape_t));

  // printf("sizeof(tape_t) = %d\n", sizeof(tape_t));
  // printf("AIMER_COMMIT_SIZE = %d\n", AIMER_COMMIT_SIZE);

  // printBytes(out_ptrs[0],sizeof(tape_t));
  // printBytes(out_ptrs[1],sizeof(tape_t));
  // printBytes(out_ptrs[2],sizeof(tape_t));
  // printBytes(out_ptrs[3],sizeof(tape_t));
  // getchar();

  // // printf("4x\n");
  // printBytes((uint8_t *)(tapes),sizeof(tape_t));
  // printBytes((uint8_t *)(tapes+1),sizeof(tape_t));
  // printBytes((uint8_t *)(tapes+2),sizeof(tape_t));
  // printBytes((uint8_t *)(tapes+3),sizeof(tape_t));
  // getchar();
}


void aim2_mpc_N2(mult_chk_N_t *mult_chk,
                const GF matrix_A[AIMER_L][AIM2_NUM_BITS_FIELD], const GF ct_GF)
{
  // pt + c = t ^ {2 ^ e - 1}
  // --> t ^ {2 ^ e} + t * c = t * pt
  // --> z = x * pt
  GF_mul_N(mult_chk->z_shares[0], (const GF *)mult_chk->x_shares[0],
           aim2_constants[0]);
  GF_transposed_matmul_add_N(mult_chk->z_shares[0],
                             (const GF *)mult_chk->x_shares[0],
                             aim2_e1_power_matrix);
  GF_transposed_matmul_add_N(mult_chk->x_shares[AIMER_L],
                             (const GF *)mult_chk->x_shares[0], matrix_A[0]);

  GF_mul_N(mult_chk->z_shares[1], (const GF *)mult_chk->x_shares[1],
           aim2_constants[1]);
  GF_transposed_matmul_add_N(mult_chk->z_shares[1],
                             (const GF *)mult_chk->x_shares[1],
                             aim2_e2_power_matrix);
  GF_transposed_matmul_add_N(mult_chk->x_shares[AIMER_L],
                             (const GF *)mult_chk->x_shares[1], matrix_A[1]);

  // x ^ {2 ^ e - 1} = pt + ct
  // --> x ^ {2 ^ e} + x * ct = x * pt
  // --> z = x * pt
  GF_sqr_N(mult_chk->z_shares[AIMER_L],
           (const GF *)mult_chk->x_shares[AIMER_L]);
  GF_sqr_N(mult_chk->z_shares[AIMER_L],
           (const GF *)mult_chk->z_shares[AIMER_L]);
  GF_sqr_N(mult_chk->z_shares[AIMER_L],
           (const GF *)mult_chk->z_shares[AIMER_L]);
  GF_mul_add_N(mult_chk->z_shares[AIMER_L],
               (const GF *)mult_chk->x_shares[AIMER_L], ct_GF);
}

void aim2_mpc_N(mult_chk_N_t *mult_chk,const GF matrix_A[AIMER_L][AIM2_NUM_BITS_FIELD], const GF ct_GF)
{
  // pt + c = t ^ {2 ^ e - 1}
  // --> t ^ {2 ^ e} + t * c = t * pt
  // --> z = x * pt
  GF_mul_N(mult_chk->z_shares[0], (const GF *)mult_chk->x_shares[0],aim2_constants[0]);
  
  BENCH_FUNCTION(GF_transposed_matmul_add_N,mult_chk->z_shares[0],(const GF *)mult_chk->x_shares[0],aim2_e1_power_matrix);
  BENCH_FUNCTION(GF_transposed_matmul_add_N_A_v1,mult_chk->z_shares[0], (const GF *)mult_chk->x_shares[0],aim2_e1_power_matrixLow,aim2_e1_power_matrixHigh);                             
  BENCH_FUNCTION(GF_transposed_matmul_add_N_A_v2,mult_chk->z_shares[0], (const GF *)mult_chk->x_shares[0],aim2_e1_power_matrixLow,aim2_e1_power_matrixHigh);
  BENCH_FUNCTION(GF_transposed_matmul_add_N_A_v3,mult_chk->z_shares[0], (const GF *)mult_chk->x_shares[0],aim2_e1_power_matrixLow,aim2_e1_power_matrixHigh);
  BENCH_FUNCTION(GF_transposed_matmul_add_N_A_v4,mult_chk->z_shares[0], (const GF *)mult_chk->x_shares[0],aim2_e1_power_matrixLowT,aim2_e1_power_matrixHighT);
  BENCH_FUNCTION(GF_transposed_matmul_add_N_A_v5,mult_chk->z_shares[0], (const GF *)mult_chk->x_shares[0],aim2_e1_power_matrixLowT,aim2_e1_power_matrixHighT);
  BENCH_FUNCTION(GF_transposed_matmul_add_N_A_v6,mult_chk->z_shares[0], (const GF *)mult_chk->x_shares[0],aim2_e1_power_matrixLowT,aim2_e1_power_matrixHighT);
  BENCH_FUNCTION(GF_transposed_matmul_add_N_A_v7,mult_chk->z_shares[0], (const GF *)mult_chk->x_shares[0],aim2_e1_power_matrixLowT,aim2_e1_power_matrixHighT);
  // BENCH_FUNCTION(GF_transposed_matmul_add_N_A,mult_chk->z_shares[0], (const GF *)mult_chk->x_shares[0],aim2_e1_power_matrixLowT,aim2_e1_power_matrixHighT);
  GF_transposed_matmul_add_N_A(mult_chk->z_shares[0], (const GF *)mult_chk->x_shares[0],aim2_e1_power_matrixLowT,aim2_e1_power_matrixHighT);  
  // /*-------------------------*/
  // vet("VET_R_AVX", (const GF *)mult_chk->z_shares[0]);
  // GF_mul_N(mult_chk->z_shares[0], (const GF *)mult_chk->x_shares[0],aim2_constants[0]);
  // vet("VET_C", mult_chk->z_shares[0]);
  // vet("VET", (const GF *)mult_chk->x_shares[0]);
  // mat("MAT", aim2_e1_power_matrix);
  // GF_transposed_matmul_add_N(mult_chk->z_shares[0],(const GF *)mult_chk->x_shares[0],aim2_e1_power_matrix);
  // vet("VET_R", mult_chk->z_shares[0]);
  // exit(1);
  // /*----------------*/

  // GF_transposed_matmul_add_N_M(mult_chk->x_shares[AIMER_L],(const GF *)mult_chk->x_shares[0], matrix_A[0]);
  // vet("VET_R", mult_chk->x_shares[AIMER_L]);
  // exit(1);
  GF_transposed_matmul_add_N(mult_chk->x_shares[AIMER_L],(const GF *)mult_chk->x_shares[0], matrix_A[0]);

  GF_mul_N(mult_chk->z_shares[1], (const GF *)mult_chk->x_shares[1],aim2_constants[1]);
  // GF_transposed_matmul_add_N(mult_chk->z_shares[1],(const GF *)mult_chk->x_shares[1],aim2_e2_power_matrix);
  GF_transposed_matmul_add_N_A(mult_chk->z_shares[1], (const GF *)mult_chk->x_shares[1],aim2_e2_power_matrixLowT,aim2_e2_power_matrixHighT);
  GF_transposed_matmul_add_N(mult_chk->x_shares[AIMER_L],(const GF *)mult_chk->x_shares[1], matrix_A[1]);

  // x ^ {2 ^ e - 1} = pt + ct
  // --> x ^ {2 ^ e} + x * ct = x * pt
  // --> z = x * pt
  GF_sqr_N(mult_chk->z_shares[AIMER_L],(const GF *)mult_chk->x_shares[AIMER_L]);
  GF_sqr_N(mult_chk->z_shares[AIMER_L],(const GF *)mult_chk->z_shares[AIMER_L]);
  GF_sqr_N(mult_chk->z_shares[AIMER_L],(const GF *)mult_chk->z_shares[AIMER_L]);
  GF_mul_add_N(mult_chk->z_shares[AIMER_L],(const GF *)mult_chk->x_shares[AIMER_L], ct_GF);
}

void aim2_mpc_N3(mult_chk_N_t *mult_chk,const GF matrix_A[AIMER_L][AIM2_NUM_BITS_FIELD], const GF ct_GF)
{
  // pt + c = t ^ {2 ^ e - 1}
  // --> t ^ {2 ^ e} + t * c = t * pt
  // --> z = x * pt
  GF_mul_N(mult_chk->z_shares[0], (const GF *)mult_chk->x_shares[0],aim2_constants[0]);
    // BENCH_FUNCTION(GF_transposed_matmul_add_N,mult_chk->z_shares[0],(const GF *)mult_chk->x_shares[0],aim2_e1_power_matrix);
  
  
  GF_transposed_matmul_add_N_A(mult_chk->z_shares[0], (const GF *)mult_chk->x_shares[0],aim2_e1_power_matrixLow,aim2_e1_power_matrixHigh);                             
  vet("VET_R_AVX", (const GF *)mult_chk->z_shares[0]);
  GF_mul_N(mult_chk->z_shares[0], (const GF *)mult_chk->x_shares[0],aim2_constants[0]);
  vet("VET_C", mult_chk->z_shares[0]);
  vet("VET", (const GF *)mult_chk->x_shares[0]);
  mat("MAT", aim2_e1_power_matrix);
  GF_transposed_matmul_add_N(mult_chk->z_shares[0],(const GF *)mult_chk->x_shares[0],aim2_e1_power_matrix);
  vet("VET_R", mult_chk->z_shares[0]);
  exit(1);

  // GF_transposed_matmul_add_N_M(mult_chk->x_shares[AIMER_L],(const GF *)mult_chk->x_shares[0], matrix_A[0]);
  // vet("VET_R", mult_chk->x_shares[AIMER_L]);
  // exit(1);
  GF_transposed_matmul_add_N(mult_chk->x_shares[AIMER_L],(const GF *)mult_chk->x_shares[0], matrix_A[0]);

  GF_mul_N(mult_chk->z_shares[1], (const GF *)mult_chk->x_shares[1],aim2_constants[1]);
  // GF_transposed_matmul_add_N(mult_chk->z_shares[1],(const GF *)mult_chk->x_shares[1],aim2_e2_power_matrix);
  GF_transposed_matmul_add_N_A(mult_chk->z_shares[1], (const GF *)mult_chk->x_shares[1],aim2_e2_power_matrixLow,aim2_e2_power_matrixHigh);
  // GF_transposed_matmul_add_N_M(mult_chk->x_shares[AIMER_L],(const GF *)mult_chk->x_shares[1], matrix_A[1]);
  GF_transposed_matmul_add_N(mult_chk->x_shares[AIMER_L],(const GF *)mult_chk->x_shares[1], matrix_A[1]);

  // x ^ {2 ^ e - 1} = pt + ct
  // --> x ^ {2 ^ e} + x * ct = x * pt
  // --> z = x * pt
  GF_sqr_N(mult_chk->z_shares[AIMER_L],(const GF *)mult_chk->x_shares[AIMER_L]);
  GF_sqr_N(mult_chk->z_shares[AIMER_L],(const GF *)mult_chk->z_shares[AIMER_L]);
  GF_sqr_N(mult_chk->z_shares[AIMER_L],(const GF *)mult_chk->z_shares[AIMER_L]);
  GF_mul_add_N(mult_chk->z_shares[AIMER_L],(const GF *)mult_chk->x_shares[AIMER_L], ct_GF);
}

// committing to the seeds and the execution views of the parties
void run_phase_1(signature_t *sign,
                 uint8_t commits[AIMER_T][AIMER_N][AIMER_COMMIT_SIZE],
                 uint8_t nodes[AIMER_T][2 * AIMER_N - 1][AIMER_SEED_SIZE],
                 mult_chk_N_t mult_chk[AIMER_T],
                 GF alpha_v_shares[AIMER_T][2][AIMER_N],
                 const uint8_t *sk, const uint8_t *m, size_t mlen)
{
  GF pt_GF = {0,}, ct_GF = {0,};
  GF_from_bytes(pt_GF, sk);
  GF_from_bytes(ct_GF, sk + AIM2_NUM_BYTES_FIELD + AIM2_IV_SIZE);

  // message pre-hashing
  hash_instance ctx;
  hash_init_prefix(&ctx, HASH_PREFIX_0);
  hash_update(&ctx, sk + AIM2_NUM_BYTES_FIELD,
              AIM2_IV_SIZE + AIM2_NUM_BYTES_FIELD);
  hash_update(&ctx, m, mlen);
  hash_final(&ctx);

  uint8_t mu[AIMER_COMMIT_SIZE];
  hash_squeeze(&ctx, mu, AIMER_COMMIT_SIZE);

  // compute first L sboxes' outputs
  GF sbox_outputs[AIMER_L];
  aim2_sbox_outputs(sbox_outputs, pt_GF);

  // derive the binary matrix and the vector from the initial vector
  GF matrix_A[AIMER_L][AIM2_NUM_BITS_FIELD];
  GF vector_b = {0,};
  generate_matrix_LU(matrix_A, vector_b, sk + AIM2_NUM_BYTES_FIELD);
  // generate_matrix_LU2(matrix_A, vector_b, sk + AIM2_NUM_BYTES_FIELD);

  // generate per-signature randomness
  uint8_t random[SECURITY_BYTES];
  randombytes(random, SECURITY_BYTES);

  // generate salt
  hash_init_prefix(&ctx, HASH_PREFIX_3);
  hash_update(&ctx, sk, AIM2_NUM_BYTES_FIELD);
  hash_update(&ctx, mu, AIMER_COMMIT_SIZE);
  hash_update(&ctx, random, SECURITY_BYTES);
  hash_final(&ctx);
  hash_squeeze(&ctx, sign->salt, AIMER_SALT_SIZE);

  // generate root seeds and expand seed trees
  for (size_t rep = 0; rep < AIMER_T; rep++)
  {
    hash_squeeze(&ctx, nodes[rep][0], AIMER_SEED_SIZE);
  }
  expand_trees(nodes, sign->salt);

  // hash_instance for h_1
  hash_init_prefix(&ctx, HASH_PREFIX_1);
  hash_update(&ctx, mu, AIMER_COMMIT_SIZE);
  hash_update(&ctx, sign->salt, AIMER_SALT_SIZE);

  hash_instance_x4 ctx_precom;
  hash_init_prefix_x4(&ctx_precom, HASH_PREFIX_5);
  
  hash_update_x4_1(&ctx_precom, sign->salt, AIMER_SALT_SIZE);

  for (size_t rep = 0; rep < AIMER_T; rep++)
  {
    // initialize adjustment values
    tape_t tapes[4];
    tape_t delta;
    memset(&delta, 0, sizeof(tape_t));

    // initialize x_star
    memset(mult_chk[rep].x_shares[AIMER_L], 0, sizeof(GF) * (AIMER_N - 1));
    GF_copy(mult_chk[rep].x_shares[AIMER_L][AIMER_N - 1], vector_b);

    for (size_t party = 0; party < AIMER_N; party += 4)
    {
      // generate execution views and commitments
      // BENCH_FUNCTION(commit_and_expand_tape_x8,tapes, commits[rep][party], sign->salt,nodes[rep][party + AIMER_N - 1], rep, party);
      // BENCH_FUNCTION(commit_and_expand_tape_x4,tapes, commits[rep][party], &ctx_precom,nodes[rep][party + AIMER_N - 1], rep, party);
      // commit_and_expand_tape_x8(tapes, commits[rep][party], sign->salt,nodes[rep][party + AIMER_N - 1], rep, party);
      commit_and_expand_tape_x4(tapes, commits[rep][party], &ctx_precom,nodes[rep][party + AIMER_N - 1], rep, party);
      hash_update(&ctx, commits[rep][party], 4 * AIMER_COMMIT_SIZE);

      for (size_t j = 0; j < 4; j++)
      {
        // compute offsets
        GF_add(delta.pt_share, delta.pt_share, tapes[j].pt_share);
        GF_add(delta.t_shares[0], delta.t_shares[0], tapes[j].t_shares[0]);
        GF_add(delta.t_shares[1], delta.t_shares[1], tapes[j].t_shares[1]);
        GF_add(delta.a_share, delta.a_share, tapes[j].a_share);
        GF_add(delta.c_share, delta.c_share, tapes[j].c_share);

        // adjust the last share and prepare the proof and h_1
        if (party + j == AIMER_N - 1)
        {
          GF_add(delta.pt_share, delta.pt_share, pt_GF);
          GF_add(delta.t_shares[0], delta.t_shares[0], sbox_outputs[0]);
          GF_add(delta.t_shares[1], delta.t_shares[1], sbox_outputs[1]);
          GF_mul_add(delta.c_share, pt_GF, delta.a_share);

          GF_to_bytes(sign->proofs[rep].delta_pt_bytes, delta.pt_share);
          GF_to_bytes(sign->proofs[rep].delta_ts_bytes[0], delta.t_shares[0]);
          GF_to_bytes(sign->proofs[rep].delta_ts_bytes[1], delta.t_shares[1]);
          GF_to_bytes(sign->proofs[rep].delta_c_bytes, delta.c_share);

          GF_add(tapes[j].pt_share, delta.pt_share, tapes[j].pt_share);
          GF_add(tapes[j].t_shares[0], delta.t_shares[0], tapes[j].t_shares[0]);
          GF_add(tapes[j].t_shares[1], delta.t_shares[1], tapes[j].t_shares[1]);
          GF_add(tapes[j].c_share, delta.c_share, tapes[j].c_share);
        }

        // run the MPC simulation and prepare the mult check inputs
        GF_copy(mult_chk[rep].pt_share[party + j], tapes[j].pt_share);
        GF_copy(mult_chk[rep].x_shares[0][party + j], tapes[j].t_shares[0]);
        GF_copy(mult_chk[rep].x_shares[1][party + j], tapes[j].t_shares[1]);
        GF_copy(alpha_v_shares[rep][0][party + j], tapes[j].a_share);
        GF_copy(alpha_v_shares[rep][1][party + j], tapes[j].c_share);
      }
    }
    // BENCH_FUNCTION(aim2_mpc_N,&mult_chk[rep], (const GF (*)[AIM2_NUM_BITS_FIELD])matrix_A,ct_GF);
    aim2_mpc_N(&mult_chk[rep], (const GF (*)[AIM2_NUM_BITS_FIELD])matrix_A,
               ct_GF);

    // NOTE: depend on the order of values in proof_t
    hash_update(&ctx, sign->proofs[rep].delta_pt_bytes,
                AIM2_NUM_BYTES_FIELD * (AIMER_L + 2));
  }

  // commit to salt, (all commitments of parties' seeds,
  // delta_pt, delta_t, delta_c) for all repetitions
  hash_final(&ctx);
  hash_squeeze(&ctx, sign->h_1, AIMER_COMMIT_SIZE);
}


// committing to the seeds and the execution views of the parties
void run_phase_1V2(signature_t *sign,
                 uint8_t commits[AIMER_T][AIMER_N][AIMER_COMMIT_SIZE],
                 uint8_t nodes[AIMER_T][2 * AIMER_N - 1][AIMER_SEED_SIZE],
                 mult_chk_N_t mult_chk[AIMER_T],
                 GF alpha_v_shares[AIMER_T][2][AIMER_N],
                 const uint8_t *sk, const uint8_t *m, size_t mlen)
{
  GF pt_GF = {0,}, ct_GF = {0,};
  GF_from_bytes(pt_GF, sk);
  GF_from_bytes(ct_GF, sk + AIM2_NUM_BYTES_FIELD + AIM2_IV_SIZE);

  // message pre-hashing
  hash_instance ctx;
  hash_init_prefix(&ctx, HASH_PREFIX_0);
  hash_update(&ctx, sk + AIM2_NUM_BYTES_FIELD,
              AIM2_IV_SIZE + AIM2_NUM_BYTES_FIELD);
  hash_update(&ctx, m, mlen);
  hash_final(&ctx);

  uint8_t mu[AIMER_COMMIT_SIZE];
  hash_squeeze(&ctx, mu, AIMER_COMMIT_SIZE);

  // compute first L sboxes' outputs
  GF sbox_outputs[AIMER_L];
  aim2_sbox_outputs(sbox_outputs, pt_GF);

  // derive the binary matrix and the vector from the initial vector
  GF matrix_A[AIMER_L][AIM2_NUM_BITS_FIELD];
  GF vector_b = {0,};
  generate_matrix_LU(matrix_A, vector_b, sk + AIM2_NUM_BYTES_FIELD);

  // generate per-signature randomness
  uint8_t random[SECURITY_BYTES];
  randombytes(random, SECURITY_BYTES);

  // generate salt
  hash_init_prefix(&ctx, HASH_PREFIX_3);
  hash_update(&ctx, sk, AIM2_NUM_BYTES_FIELD);
  hash_update(&ctx, mu, AIMER_COMMIT_SIZE);
  hash_update(&ctx, random, SECURITY_BYTES);
  hash_final(&ctx);
  hash_squeeze(&ctx, sign->salt, AIMER_SALT_SIZE);

  // generate root seeds and expand seed trees
  for (size_t rep = 0; rep < AIMER_T; rep++)
  {
    hash_squeeze(&ctx, nodes[rep][0], AIMER_SEED_SIZE);
  }
  expand_trees(nodes, sign->salt);

  // hash_instance for h_1
  hash_init_prefix(&ctx, HASH_PREFIX_1);
  hash_update(&ctx, mu, AIMER_COMMIT_SIZE);
  hash_update(&ctx, sign->salt, AIMER_SALT_SIZE);

  hash_instance_x4 ctx_precom;
  hash_init_prefix_x4(&ctx_precom, HASH_PREFIX_5);
  
  hash_update_x4_1(&ctx_precom, sign->salt, AIMER_SALT_SIZE);

  for (size_t rep = 0; rep < AIMER_T; rep++)
  {
    // initialize adjustment values
    tape_t tapes[8];
    tape_t delta;
    memset(&delta, 0, sizeof(tape_t));

    // initialize x_star
    memset(mult_chk[rep].x_shares[AIMER_L], 0, sizeof(GF) * (AIMER_N - 1));
    GF_copy(mult_chk[rep].x_shares[AIMER_L][AIMER_N - 1], vector_b);

    for (size_t party = 0; party < AIMER_N; party += 8)
    {
      // generate execution views and commitments
      // BENCH_FUNCTION(commit_and_expand_tape_x8,tapes, commits[rep][party], sign->salt,nodes[rep][party + AIMER_N - 1], rep, party);
      commit_and_expand_tape_x8(tapes, commits[rep][party], sign->salt,nodes[rep][party + AIMER_N - 1], rep, party);
      // commit_and_expand_tape_x4(tapes, commits[rep][party], &ctx_precom,nodes[rep][party + AIMER_N - 1], rep, party);
      hash_update(&ctx, commits[rep][party], 8 * AIMER_COMMIT_SIZE);

      for (size_t j = 0; j < 8; j++)
      {
        // compute offsets
        GF_add(delta.pt_share, delta.pt_share, tapes[j].pt_share);
        GF_add(delta.t_shares[0], delta.t_shares[0], tapes[j].t_shares[0]);
        GF_add(delta.t_shares[1], delta.t_shares[1], tapes[j].t_shares[1]);
        GF_add(delta.a_share, delta.a_share, tapes[j].a_share);
        GF_add(delta.c_share, delta.c_share, tapes[j].c_share);

        // adjust the last share and prepare the proof and h_1
        if (party + j == AIMER_N - 1)
        {
          GF_add(delta.pt_share, delta.pt_share, pt_GF);
          GF_add(delta.t_shares[0], delta.t_shares[0], sbox_outputs[0]);
          GF_add(delta.t_shares[1], delta.t_shares[1], sbox_outputs[1]);
          GF_mul_add(delta.c_share, pt_GF, delta.a_share);

          GF_to_bytes(sign->proofs[rep].delta_pt_bytes, delta.pt_share);
          GF_to_bytes(sign->proofs[rep].delta_ts_bytes[0], delta.t_shares[0]);
          GF_to_bytes(sign->proofs[rep].delta_ts_bytes[1], delta.t_shares[1]);
          GF_to_bytes(sign->proofs[rep].delta_c_bytes, delta.c_share);

          GF_add(tapes[j].pt_share, delta.pt_share, tapes[j].pt_share);
          GF_add(tapes[j].t_shares[0], delta.t_shares[0], tapes[j].t_shares[0]);
          GF_add(tapes[j].t_shares[1], delta.t_shares[1], tapes[j].t_shares[1]);
          GF_add(tapes[j].c_share, delta.c_share, tapes[j].c_share);
        }

        // run the MPC simulation and prepare the mult check inputs
        GF_copy(mult_chk[rep].pt_share[party + j], tapes[j].pt_share);
        GF_copy(mult_chk[rep].x_shares[0][party + j], tapes[j].t_shares[0]);
        GF_copy(mult_chk[rep].x_shares[1][party + j], tapes[j].t_shares[1]);
        GF_copy(alpha_v_shares[rep][0][party + j], tapes[j].a_share);
        GF_copy(alpha_v_shares[rep][1][party + j], tapes[j].c_share);
      }
    }
    aim2_mpc_N(&mult_chk[rep], (const GF (*)[AIM2_NUM_BITS_FIELD])matrix_A,ct_GF);

    // NOTE: depend on the order of values in proof_t
    hash_update(&ctx, sign->proofs[rep].delta_pt_bytes,AIM2_NUM_BYTES_FIELD * (AIMER_L + 2));
  }

  // commit to salt, (all commitments of parties' seeds,
  // delta_pt, delta_t, delta_c) for all repetitions
  hash_final(&ctx);
  hash_squeeze(&ctx, sign->h_1, AIMER_COMMIT_SIZE);
}

void run_phase_2_and_3(signature_t *sign,
                       GF alpha_v_shares[AIMER_T][2][AIMER_N],
                       const mult_chk_N_t mult_chk[AIMER_T])
{
  hash_instance ctx_e;
  hash_init(&ctx_e);
  hash_update(&ctx_e, sign->h_1, AIMER_COMMIT_SIZE);
  hash_final(&ctx_e);

  hash_instance ctx;
  hash_init_prefix(&ctx, HASH_PREFIX_2);
  hash_update(&ctx, sign->h_1, AIMER_COMMIT_SIZE);
  hash_update(&ctx, sign->salt, AIMER_SALT_SIZE);

  for (size_t rep = 0; rep < AIMER_T; rep++)
  {
    GF alpha = {0,};
    GF epsilons[AIMER_L + 1];
    hash_squeeze(&ctx_e, (uint8_t *)epsilons, sizeof(epsilons));

    GF alpha_v_shares_hi[2][AIMER_N];
    memset(alpha_v_shares_hi, 0, sizeof(alpha_v_shares_hi));

    POLY_mul_add_N(alpha_v_shares[rep][0], alpha_v_shares_hi[0],
                   mult_chk[rep].x_shares[0], epsilons[0]);
    POLY_mul_add_N(alpha_v_shares[rep][1], alpha_v_shares_hi[1],
                   mult_chk[rep].z_shares[0], epsilons[0]);

    POLY_mul_add_N(alpha_v_shares[rep][0], alpha_v_shares_hi[0],
                   mult_chk[rep].x_shares[1], epsilons[1]);
    POLY_mul_add_N(alpha_v_shares[rep][1], alpha_v_shares_hi[1],
                   mult_chk[rep].z_shares[1], epsilons[1]);

    POLY_mul_add_N(alpha_v_shares[rep][0], alpha_v_shares_hi[0],
                   mult_chk[rep].x_shares[2], epsilons[2]);
    POLY_mul_add_N(alpha_v_shares[rep][1], alpha_v_shares_hi[1],
                   mult_chk[rep].z_shares[2], epsilons[2]);

    POLY_red_N(alpha_v_shares[rep][0], (const GF *)alpha_v_shares_hi[0]);
    POLY_red_N(alpha_v_shares[rep][1], (const GF *)alpha_v_shares_hi[1]);

    for (size_t party = 0; party < AIMER_N; party++)
    {
      GF_add(alpha, alpha, alpha_v_shares[rep][0][party]);
    }

    // alpha is opened, so we can finish calculating v_share
    GF_mul_add_N(alpha_v_shares[rep][1], mult_chk[rep].pt_share, alpha);
    hash_update(&ctx, (const uint8_t *)alpha_v_shares[rep],
                AIM2_NUM_BYTES_FIELD * 2 * AIMER_N);
  }

  hash_final(&ctx);
  hash_squeeze(&ctx, sign->h_2, AIMER_COMMIT_SIZE);
}

////////////////////////////////////////////////////////////////////////////////
int crypto_sign_keypair(uint8_t *pk, uint8_t *sk)
{
  if (!pk || !sk)
  {
    return -1;
  }

  randombytes(sk, AIM2_NUM_BYTES_FIELD);
  randombytes(pk, AIM2_IV_SIZE);

  aim2(pk + AIM2_IV_SIZE, sk, pk);
  memcpy(sk + AIM2_NUM_BYTES_FIELD, pk, AIM2_IV_SIZE + AIM2_NUM_BYTES_FIELD);

  return 0;
}

int crypto_sign_signature(uint8_t *sig, size_t *siglen,
        const uint8_t *m, size_t mlen,
        const uint8_t *sk)
{
  hash_instance ctx;
  signature_t *sign = (signature_t *)sig;

  //////////////////////////////////////////////////////////////////////////
  // Phase 1: Committing to the seeds and the execution views of parties. //
  //////////////////////////////////////////////////////////////////////////

  // nodes for seed trees
  uint8_t nodes[AIMER_T][2 * AIMER_N - 1][AIMER_SEED_SIZE];

  // commitments for seeds
  uint8_t commits[AIMER_T][AIMER_N][AIMER_COMMIT_SIZE];

  // multiplication check inputs
  mult_chk_N_t mult_chk[AIMER_T];

  // multiplication check outputs
  GF alpha_v_shares[AIMER_T][2][AIMER_N];

  #ifdef benchPhases
                  BENCH_FUNCTION(run_phase_1,sign, commits, (uint8_t (*)[2 * AIMER_N - 1][AIMER_SEED_SIZE])nodes, mult_chk, alpha_v_shares, sk, m, mlen);  
  #endif

  // BENCH_FUNCTION(run_phase_1V2,sign, commits,(uint8_t (*)[2 * AIMER_N - 1][AIMER_SEED_SIZE])nodes, mult_chk,alpha_v_shares, sk, m, mlen);
  // BENCH_FUNCTION(run_phase_1,sign, commits, (uint8_t (*)[2 * AIMER_N - 1][AIMER_SEED_SIZE])nodes, mult_chk, alpha_v_shares, sk, m, mlen);
  // run_phase_1V2(sign, commits,(uint8_t (*)[2 * AIMER_N - 1][AIMER_SEED_SIZE])nodes, mult_chk,alpha_v_shares, sk, m, mlen);
  run_phase_1(sign, commits,(uint8_t (*)[2 * AIMER_N - 1][AIMER_SEED_SIZE])nodes, mult_chk,alpha_v_shares, sk, m, mlen);

  /////////////////////////////////////////////////////////////////
  // Phase 2, 3: Challenging and committing to the simulation of //
  //             the multiplication checking protocol.           //
  /////////////////////////////////////////////////////////////////

  // compute the commitment of phase 3
  #ifdef benchPhases
                  BENCH_FUNCTION(run_phase_2_and_3,sign, alpha_v_shares, mult_chk);
  #endif            
  run_phase_2_and_3(sign, alpha_v_shares, mult_chk);

  //////////////////////////////////////////////////////
  // Phase 4: Challenging views of the MPC protocols. //
  //////////////////////////////////////////////////////

  hash_init(&ctx);
  hash_update(&ctx, sign->h_2, AIMER_COMMIT_SIZE);
  hash_final(&ctx);

  uint8_t indices[AIMER_T]; // AIMER_N <= 256
  hash_squeeze(&ctx, indices, AIMER_T);
  for (size_t rep = 0; rep < AIMER_T; rep++)
  {
    indices[rep] &= (1 << AIMER_LOGN) - 1;
  }

  //////////////////////////////////////////////////////
  // Phase 5: Opening the views of the MPC protocols. //
  //////////////////////////////////////////////////////

  crypto_declassify(indices, sizeof(indices));
  for (size_t rep = 0; rep < AIMER_T; rep++)
  {
    size_t i_bar = indices[rep];
    reveal_all_but(sign->proofs[rep].reveal_path,
                   (const uint8_t (*)[AIMER_SEED_SIZE])nodes[rep], i_bar);
    memcpy(sign->proofs[rep].missing_commitment, commits[rep][i_bar],
           AIMER_COMMIT_SIZE);
    GF_to_bytes(sign->proofs[rep].missing_alpha_share_bytes,
                alpha_v_shares[rep][0][i_bar]);
  }
  *siglen = CRYPTO_BYTES;

  return 0;
}

int crypto_sign(uint8_t *sm, size_t *smlen,
        const uint8_t *m, size_t mlen,
        const uint8_t *sk)
{
  crypto_sign_signature(sm + mlen, smlen, m, mlen, sk);

  memcpy(sm, m, mlen);
  *smlen += mlen;

  return 0;
}

int crypto_sign_verify(const uint8_t *sig, size_t siglen,
        const uint8_t *m, size_t mlen,
        const uint8_t *pk)
{
  if (siglen != CRYPTO_BYTES)
  {
    return -1;
  }

  const signature_t *sign = (const signature_t *)sig;

  GF ct_GF = {0,};
  GF_from_bytes(ct_GF, pk + AIM2_IV_SIZE);

  // derive the binary matrix and the vector from the initial vector
  GF matrix_A[AIMER_L][AIM2_NUM_BITS_FIELD];
  GF vector_b = {0,};
  generate_matrix_LU(matrix_A, vector_b, pk);

  hash_instance ctx_e, ctx_h1, ctx_h2;

  // indices = Expand(h_2)
  hash_init(&ctx_e);
  hash_update(&ctx_e, sign->h_2, AIMER_COMMIT_SIZE);
  hash_final(&ctx_e);

  uint8_t indices[AIMER_T]; // AIMER_N <= 256
  hash_squeeze(&ctx_e, indices, AIMER_T);
  for (size_t rep = 0; rep < AIMER_T; rep++)
  {
    indices[rep] &= (1 << AIMER_LOGN) - 1;
  }

  // epsilons = Expand(h_1)
  hash_init(&ctx_e);
  hash_update(&ctx_e, sign->h_1, AIMER_COMMIT_SIZE);
  hash_final(&ctx_e);

  // message pre-hashing
  uint8_t mu[AIMER_COMMIT_SIZE];
  hash_init_prefix(&ctx_h1, HASH_PREFIX_0);
  hash_update(&ctx_h1, pk, AIM2_IV_SIZE + AIM2_NUM_BYTES_FIELD);
  hash_update(&ctx_h1, m, mlen);
  hash_final(&ctx_h1);
  hash_squeeze(&ctx_h1, mu, AIMER_COMMIT_SIZE);

  // ready for computing h_1' and h_2'
  hash_init_prefix(&ctx_h1, HASH_PREFIX_1);
  hash_update(&ctx_h1, mu, AIMER_COMMIT_SIZE);
  hash_update(&ctx_h1, sign->salt, AIMER_SALT_SIZE);

  hash_init_prefix(&ctx_h2, HASH_PREFIX_2);
  hash_update(&ctx_h2, sign->h_1, AIMER_COMMIT_SIZE);
  hash_update(&ctx_h2, sign->salt, AIMER_SALT_SIZE);

  hash_instance_x4 ctx_precom;
  hash_init_prefix_x4(&ctx_precom, HASH_PREFIX_5);
  hash_update_x4_1(&ctx_precom, sign->salt, AIMER_SALT_SIZE);

  for (size_t rep = 0; rep < AIMER_T; rep++)
  {
    size_t i_bar = indices[rep];
    uint8_t nodes[2 * AIMER_N - 2][AIMER_SEED_SIZE];
    memset(nodes, 0, sizeof(nodes));

    reconstruct_tree(nodes, sign->salt, sign->proofs[rep].reveal_path,
                     rep, i_bar);

    mult_chk_N_t mult_chk;
    memset(&mult_chk, 0, sizeof(mult_chk_N_t));
    GF_copy(mult_chk.x_shares[AIMER_L][AIMER_N - 1], vector_b);

    GF epsilons[AIMER_L + 1];
    hash_squeeze(&ctx_e, (uint8_t *)epsilons, sizeof(epsilons));

    GF alpha_v = {0,};
    GF alpha_v_shares[2][AIMER_N], alpha_v_shares_hi[2][AIMER_N];
    memset(alpha_v_shares_hi, 0, sizeof(alpha_v_shares_hi));

    for (size_t party = 0; party < AIMER_N; party += 4)
    {
      tape_t tapes[4];
      uint8_t commits[4][AIMER_COMMIT_SIZE];
      commit_and_expand_tape_x4(tapes, commits[0], &ctx_precom,
                                nodes[party + AIMER_N - 2], rep, party);

      if (party / 4 == i_bar / 4)
      {
        memcpy(commits[i_bar % 4], sign->proofs[rep].missing_commitment,
               AIMER_COMMIT_SIZE);
      }
      hash_update(&ctx_h1, commits[0], 4 * AIMER_COMMIT_SIZE);

      if (party == AIMER_N - 4)
      {
        GF temp = {0,};

        GF_from_bytes(temp, sign->proofs[rep].delta_pt_bytes);
        GF_add(tapes[3].pt_share, tapes[3].pt_share, temp);

        GF_from_bytes(temp, sign->proofs[rep].delta_ts_bytes[0]);
        GF_add(tapes[3].t_shares[0], tapes[3].t_shares[0], temp);
        GF_from_bytes(temp, sign->proofs[rep].delta_ts_bytes[1]);
        GF_add(tapes[3].t_shares[1], tapes[3].t_shares[1], temp);

        GF_from_bytes(temp, sign->proofs[rep].delta_c_bytes);
        GF_add(tapes[3].c_share, tapes[3].c_share, temp);
      }

      for (size_t j = 0; j < 4; j++)
      {
        GF_copy(mult_chk.pt_share[party + j], tapes[j].pt_share);
        GF_copy(mult_chk.x_shares[0][party + j], tapes[j].t_shares[0]);
        GF_copy(mult_chk.x_shares[1][party + j], tapes[j].t_shares[1]);
        GF_copy(alpha_v_shares[0][party + j], tapes[j].a_share);
        GF_copy(alpha_v_shares[1][party + j], tapes[j].c_share);
      }
    }

    aim2_mpc_N(&mult_chk, (const GF (*)[AIM2_NUM_BITS_FIELD])matrix_A, ct_GF);
    hash_update(&ctx_h1, sign->proofs[rep].delta_pt_bytes,
                AIM2_NUM_BYTES_FIELD * (AIMER_L + 2));

    POLY_mul_add_N(alpha_v_shares[0], alpha_v_shares_hi[0],
                   (const GF *)mult_chk.x_shares[0], epsilons[0]);
    POLY_mul_add_N(alpha_v_shares[1], alpha_v_shares_hi[1],
                   (const GF *)mult_chk.z_shares[0], epsilons[0]);

    POLY_mul_add_N(alpha_v_shares[0], alpha_v_shares_hi[0],
                   (const GF *)mult_chk.x_shares[1], epsilons[1]);
    POLY_mul_add_N(alpha_v_shares[1], alpha_v_shares_hi[1],
                   (const GF *)mult_chk.z_shares[1], epsilons[1]);

    POLY_mul_add_N(alpha_v_shares[0], alpha_v_shares_hi[0],
                   (const GF *)mult_chk.x_shares[2], epsilons[2]);
    POLY_mul_add_N(alpha_v_shares[1], alpha_v_shares_hi[1],
                   (const GF *)mult_chk.z_shares[2], epsilons[2]);

    POLY_red_N(alpha_v_shares[0], (const GF *)alpha_v_shares_hi[0]);
    POLY_red_N(alpha_v_shares[1], (const GF *)alpha_v_shares_hi[1]);

    // open alpha
    GF_from_bytes(alpha_v_shares[0][i_bar],
                  sign->proofs[rep].missing_alpha_share_bytes);
    GF_set0(alpha_v);
    for (size_t party = 0; party < AIMER_N; party++)
    {
      GF_add(alpha_v, alpha_v, alpha_v_shares[0][party]);
    }

    // alpha is opened, so we can finish calculating v_share
    GF_mul_add_N(alpha_v_shares[1], (const GF *)mult_chk.pt_share, alpha_v);
    GF_set0(alpha_v);
    for (size_t party = 0; party < AIMER_N; party++)
    {
      GF_add(alpha_v, alpha_v, alpha_v_shares[1][party]);
    }
    GF_add(alpha_v_shares[1][i_bar], alpha_v, alpha_v_shares[1][i_bar]);
    hash_update(&ctx_h2, (const uint8_t *)alpha_v_shares,
                AIM2_NUM_BYTES_FIELD * 2 * AIMER_N);
  }

  uint8_t h_1_prime[AIMER_COMMIT_SIZE];
  hash_final(&ctx_h1);
  hash_squeeze(&ctx_h1, h_1_prime, AIMER_COMMIT_SIZE);

  uint8_t h_2_prime[AIMER_COMMIT_SIZE];
  hash_final(&ctx_h2);
  hash_squeeze(&ctx_h2, h_2_prime, AIMER_COMMIT_SIZE);

  if (memcmp(h_1_prime, sign->h_1, AIMER_COMMIT_SIZE) != 0 ||
      memcmp(h_2_prime, sign->h_2, AIMER_COMMIT_SIZE) != 0)
  {
    return -1;
  }

  return 0;
}

int crypto_sign_open(uint8_t *m, size_t *mlen,
        const uint8_t *sm, size_t smlen,
        const uint8_t *pk)
{
  if (smlen < CRYPTO_BYTES)
  {
    return -1;
  }

  const size_t message_len = smlen - CRYPTO_BYTES;
  const uint8_t *message = sm;
  const uint8_t *signature = sm + message_len;

  if (crypto_sign_verify(signature, CRYPTO_BYTES, message, message_len, pk))
  {
    return -1;
  }

  memcpy(m, message, message_len);
  *mlen = message_len;

  return 0;
}
