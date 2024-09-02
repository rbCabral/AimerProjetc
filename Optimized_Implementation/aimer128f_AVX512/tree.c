// SPDX-License-Identifier: MIT

#include "tree.h"
#include "hash.h"
#include <stddef.h>
#include <stdint.h>
#include <string.h>

//  Example of tree for [N = 8]
//  x
//  d = 0: 1
//  d = 1: 2         3
//  d = 2: 4   5     6     7
//  d = 3: 8 9 10 11 12 13 14 15



void expand_trees_AVX512(uint8_t nodes[AIMER_T][2 * AIMER_N - 1][AIMER_SEED_SIZE],const uint8_t salt[AIMER_SALT_SIZE])
{
  uint8_t bufs[4][AIMER_SEED_SIZE + 2];
  const uint8_t *in_ptrs[4] = {bufs[0], bufs[1], bufs[2], bufs[3]};
  uint8_t *out_ptrs[4];

  hash_instance_x4 ctx, ctx_;
  hash_init_prefix_x4(&ctx_, HASH_PREFIX_4);
  hash_update_x4_1(&ctx_, salt, AIMER_SALT_SIZE);
  // uint8_t out[64];  

  __m512i s[25],aux0,aux1;
  uint8_t bufs2[4][64];
  uint8_t init[17];
  const __m512i m0 = _mm512_set_epi64(0x00,0xFFFFFFFF,0xFFFFFFFFFFFFFFFF,0xFFFFFFFFFFFFFFFF,0x00,0xFFFFFFFF,0xFFFFFFFFFFFFFFFF,0xFFFFFFFFFFFFFFFF);
  
  init[0] = 0x04;
  memcpy(init+1, salt,16);

  aux0 = _mm512_set1_epi64(((uint64_t*)init)[0]);
  aux1 = _mm512_set1_epi64(((uint64_t*)init)[1]);

  size_t rep, index, depth;
  for (rep = 0; rep + 8 <= AIMER_T; rep += 8)
  {
    for (index = 1; index < AIMER_N; index++)
    {
      for(int i=0;i<4;i++){
        bufs2[i][0]  = init[16];//last byte of the salt  
        bufs2[i][32] = init[16];  
        bufs2[i][1]  = (uint8_t)(rep+i);//repetition
        bufs2[i][33] = (uint8_t)(rep+ i + 4);
        bufs2[i][2]  = (uint8_t)(index);//party
        bufs2[i][34] = (uint8_t)(index);
        memcpy(&bufs2[i][3], nodes[rep + i][index - 1], AIMER_SEED_SIZE);//seed
        memcpy(&bufs2[i][35], nodes[rep + i+ 4][index - 1] , AIMER_SEED_SIZE);
        bufs2[i][19] = 0x1F;//shake128
        bufs2[i][51] = 0x1F;  
      }
      s[0] = aux0;
      s[1] = aux1;
      s[2] = _mm512_and_si512(_mm512_loadu_si512((__m512i*)bufs2[0]),m0);
      s[3] = _mm512_and_si512(_mm512_loadu_si512((__m512i*)bufs2[1]),m0);
      s[4] = _mm512_and_si512(_mm512_loadu_si512((__m512i*)bufs2[2]),m0);
      s[5] = _mm512_and_si512(_mm512_loadu_si512((__m512i*)bufs2[3]),m0);

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
        // _mm512_storeu_si512((__m512i*)out,s[i]);
        // printBytes(out,64);
        // }
    
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

      _mm256_storeu_si256((__m256i*)nodes[rep + 0][2 * index - 1],_mm512_extracti64x4_epi64(s[19],0));
      _mm256_storeu_si256((__m256i*)nodes[rep + 1][2 * index - 1],_mm512_extracti64x4_epi64(s[20],0));
      _mm256_storeu_si256((__m256i*)nodes[rep + 2][2 * index - 1],_mm512_extracti64x4_epi64(s[19],1));
      _mm256_storeu_si256((__m256i*)nodes[rep + 3][2 * index - 1],_mm512_extracti64x4_epi64(s[20],1));
      _mm256_storeu_si256((__m256i*)nodes[rep + 4][2 * index - 1],_mm512_extracti64x4_epi64(s[21],0));
      _mm256_storeu_si256((__m256i*)nodes[rep + 5][2 * index - 1],_mm512_extracti64x4_epi64(s[22],0));
      _mm256_storeu_si256((__m256i*)nodes[rep + 6][2 * index - 1],_mm512_extracti64x4_epi64(s[21],1));
      _mm256_storeu_si256((__m256i*)nodes[rep + 7][2 * index - 1],_mm512_extracti64x4_epi64(s[22],1));
      

  // printBytes((uint8_t *)nodes[rep + 0][2 * index - 1],32);
  // printBytes((uint8_t *)nodes[rep + 1][2 * index - 1],32);
  // printBytes((uint8_t *)nodes[rep + 2][2 * index - 1],32);
  // printBytes((uint8_t *)nodes[rep + 3][2 * index - 1],32);
  // printBytes((uint8_t *)nodes[rep + 4][2 * index - 1],32);
  // printBytes((uint8_t *)nodes[rep + 5][2 * index - 1],32);
  // printBytes((uint8_t *)nodes[rep + 6][2 * index - 1],32);
  // printBytes((uint8_t *)nodes[rep + 7][2 * index - 1],32);
  // getchar();

      // hash_update_x4(&ctx, in_ptrs, AIMER_SEED_SIZE + 2);
      // hash_final_x4_2(&ctx);

      // out_ptrs[0] = nodes[rep + 0][2 * index - 1];
      // out_ptrs[1] = nodes[rep + 1][2 * index - 1];
      // out_ptrs[2] = nodes[rep + 2][2 * index - 1];
      // out_ptrs[3] = nodes[rep + 3][2 * index - 1];

      // hash_squeeze_x4(&ctx, out_ptrs, AIMER_SEED_SIZE << 1);
    }
  }

  hash_instance ctx1, ctx1_;
  hash_init_prefix(&ctx1_, HASH_PREFIX_4);
  hash_update(&ctx1_, salt, AIMER_SALT_SIZE);

  for (; rep < AIMER_T; rep++)
  {
    for (depth = 0; depth < 2; depth++)
    {
      for (index = (1U << depth); index < (2U << depth); index++) 
      {
        memcpy(&ctx1, &ctx1_, sizeof(hash_instance));
        bufs[0][0] = (uint8_t)(rep);
        bufs[0][1] = (uint8_t)(index);
        memcpy(bufs[0] + 2, nodes[rep][index - 1], AIMER_SEED_SIZE);

        hash_update(&ctx1, bufs[0], AIMER_SEED_SIZE + 2);
        hash_final(&ctx1);
        hash_squeeze(&ctx1, nodes[rep][2 * index - 1], AIMER_SEED_SIZE << 1);
      }
    }
    // depth = 2;
    for (; depth < AIMER_LOGN; depth++)
    {
      int p=0;
      for (index = (1U << depth); index < (2U << depth); index += 4)
      {
        p++;
        memcpy(&ctx, &ctx_, sizeof(hash_instance_x4));

        bufs[0][0] = (uint8_t)(rep);
        bufs[1][0] = (uint8_t)(rep);
        bufs[2][0] = (uint8_t)(rep);
        bufs[3][0] = (uint8_t)(rep);

        bufs[0][1] = (uint8_t)(index + 0);
        bufs[1][1] = (uint8_t)(index + 1);
        bufs[2][1] = (uint8_t)(index + 2);
        bufs[3][1] = (uint8_t)(index + 3);

        memcpy(bufs[0] + 2, nodes[rep][index - 1], AIMER_SEED_SIZE);
        memcpy(bufs[1] + 2, nodes[rep][index + 0], AIMER_SEED_SIZE);
        memcpy(bufs[2] + 2, nodes[rep][index + 1], AIMER_SEED_SIZE);
        memcpy(bufs[3] + 2, nodes[rep][index + 2], AIMER_SEED_SIZE);

        hash_update_x4(&ctx, in_ptrs, AIMER_SEED_SIZE + 2);
        hash_final_x4(&ctx);

        out_ptrs[0] = nodes[rep][2 * index - 1];
        out_ptrs[1] = nodes[rep][2 * index + 1];
        out_ptrs[2] = nodes[rep][2 * index + 3];
        out_ptrs[3] = nodes[rep][2 * index + 5];

        hash_squeeze_x4(&ctx, out_ptrs, AIMER_SEED_SIZE << 1);
      }
      // printf("p = %d\n",p);
    }
    // getchar();
  }
}

void expand_trees(uint8_t nodes[AIMER_T][2 * AIMER_N - 1][AIMER_SEED_SIZE],
                  const uint8_t salt[AIMER_SALT_SIZE])
{
  uint8_t bufs[4][AIMER_SEED_SIZE + 2];
  const uint8_t *in_ptrs[4] = {bufs[0], bufs[1], bufs[2], bufs[3]};
  uint8_t *out_ptrs[4];

  hash_instance_x4 ctx, ctx_;
  hash_init_prefix_x4(&ctx_, HASH_PREFIX_4);
  hash_update_x4_1(&ctx_, salt, AIMER_SALT_SIZE);

  size_t rep, index, depth;
  for (rep = 0; rep + 4 <= AIMER_T; rep += 4)
  {
    for (index = 1; index < AIMER_N; index++)
    {
      memcpy(&ctx, &ctx_, sizeof(hash_instance_x4));

      bufs[0][0] = (uint8_t)(rep + 0);
      bufs[1][0] = (uint8_t)(rep + 1);
      bufs[2][0] = (uint8_t)(rep + 2);
      bufs[3][0] = (uint8_t)(rep + 3);

      bufs[0][1] = (uint8_t)(index);
      bufs[1][1] = (uint8_t)(index);
      bufs[2][1] = (uint8_t)(index);
      bufs[3][1] = (uint8_t)(index);

      memcpy(bufs[0] + 2, nodes[rep + 0][index - 1], AIMER_SEED_SIZE);
      memcpy(bufs[1] + 2, nodes[rep + 1][index - 1], AIMER_SEED_SIZE);
      memcpy(bufs[2] + 2, nodes[rep + 2][index - 1], AIMER_SEED_SIZE);
      memcpy(bufs[3] + 2, nodes[rep + 3][index - 1], AIMER_SEED_SIZE);

      hash_update_x4(&ctx, in_ptrs, AIMER_SEED_SIZE + 2);
      hash_final_x4(&ctx);

      out_ptrs[0] = nodes[rep + 0][2 * index - 1];
      out_ptrs[1] = nodes[rep + 1][2 * index - 1];
      out_ptrs[2] = nodes[rep + 2][2 * index - 1];
      out_ptrs[3] = nodes[rep + 3][2 * index - 1];

      hash_squeeze_x4(&ctx, out_ptrs, AIMER_SEED_SIZE << 1);
    }
  }

  hash_instance ctx1, ctx1_;
  hash_init_prefix(&ctx1_, HASH_PREFIX_4);
  hash_update(&ctx1_, salt, AIMER_SALT_SIZE);

  for (; rep < AIMER_T; rep++)
  {
    for (depth = 0; depth < 2; depth++)
    {
      for (index = (1U << depth); index < (2U << depth); index++) 
      {
        memcpy(&ctx1, &ctx1_, sizeof(hash_instance));
        bufs[0][0] = (uint8_t)(rep);
        bufs[0][1] = (uint8_t)(index);
        memcpy(bufs[0] + 2, nodes[rep][index - 1], AIMER_SEED_SIZE);

        hash_update(&ctx1, bufs[0], AIMER_SEED_SIZE + 2);
        hash_final(&ctx1);
        hash_squeeze(&ctx1, nodes[rep][2 * index - 1], AIMER_SEED_SIZE << 1);
      }
    }
    // depth = 2;
    for (; depth < AIMER_LOGN; depth++)
    {
      for (index = (1U << depth); index < (2U << depth); index += 4)
      {
        memcpy(&ctx, &ctx_, sizeof(hash_instance_x4));

        bufs[0][0] = (uint8_t)(rep);
        bufs[1][0] = (uint8_t)(rep);
        bufs[2][0] = (uint8_t)(rep);
        bufs[3][0] = (uint8_t)(rep);

        bufs[0][1] = (uint8_t)(index + 0);
        bufs[1][1] = (uint8_t)(index + 1);
        bufs[2][1] = (uint8_t)(index + 2);
        bufs[3][1] = (uint8_t)(index + 3);

        memcpy(bufs[0] + 2, nodes[rep][index - 1], AIMER_SEED_SIZE);
        memcpy(bufs[1] + 2, nodes[rep][index + 0], AIMER_SEED_SIZE);
        memcpy(bufs[2] + 2, nodes[rep][index + 1], AIMER_SEED_SIZE);
        memcpy(bufs[3] + 2, nodes[rep][index + 2], AIMER_SEED_SIZE);

        hash_update_x4(&ctx, in_ptrs, AIMER_SEED_SIZE + 2);
        hash_final_x4(&ctx);

        out_ptrs[0] = nodes[rep][2 * index - 1];
        out_ptrs[1] = nodes[rep][2 * index + 1];
        out_ptrs[2] = nodes[rep][2 * index + 3];
        out_ptrs[3] = nodes[rep][2 * index + 5];

        hash_squeeze_x4(&ctx, out_ptrs, AIMER_SEED_SIZE << 1);
      }
    }
  }
}

void reveal_all_but(uint8_t reveal_path[AIMER_LOGN][AIMER_SEED_SIZE],
                    const uint8_t nodes[2 * AIMER_N - 1][AIMER_SEED_SIZE],
                    size_t cover_index)
{
  size_t index = cover_index + AIMER_N;
  for (size_t depth = 0; depth < AIMER_LOGN; depth++)
  {
    // index ^ 1 is sibling index
    memcpy(reveal_path[depth], nodes[(index ^ 1) - 1], AIMER_SEED_SIZE);

    // go to parent node
    index >>= 1;
  }
}

void reconstruct_tree(uint8_t nodes[2 * AIMER_N - 2][AIMER_SEED_SIZE],
                      const uint8_t salt[AIMER_SALT_SIZE],
                      const uint8_t reveal_path[AIMER_LOGN][AIMER_SEED_SIZE],
                      size_t rep_index,
                      size_t cover_index)
{
  size_t index, depth, path;
  uint8_t bufs[4][AIMER_SEED_SIZE + 2];
  const uint8_t *in_ptrs[4] = {bufs[0], bufs[1], bufs[2], bufs[3]};
  uint8_t *out_ptrs[4];

  hash_instance_x4 ctx, ctx_;
  hash_init_prefix_x4(&ctx_, HASH_PREFIX_4);
  hash_update_x4_1(&ctx_, salt, AIMER_SALT_SIZE);

  // depth = 1
  hash_instance ctx1;
  uint8_t buffer[AIMER_SALT_SIZE + AIMER_SEED_SIZE + 3];

  hash_init(&ctx1);
  path = ((cover_index + AIMER_N) >> (AIMER_LOGN - 1)) ^ 1;

  buffer[0] = HASH_PREFIX_4;
  memcpy(buffer + 1, salt, AIMER_SALT_SIZE);
  buffer[AIMER_SALT_SIZE + 1] = rep_index;
  buffer[AIMER_SALT_SIZE + 2] = path;
  memcpy(buffer + AIMER_SALT_SIZE + 3, reveal_path[AIMER_LOGN - 1],
         AIMER_SEED_SIZE);

  hash_update(&ctx1, buffer, AIMER_SALT_SIZE+ AIMER_SEED_SIZE + 3);
  hash_final(&ctx1);
  hash_squeeze(&ctx1, nodes[(path << 1) - 2], 2 * AIMER_SEED_SIZE);

  for (depth = 2; depth < AIMER_LOGN; depth++)
  {
    path = ((cover_index + AIMER_N) >> (AIMER_LOGN - depth)) ^ 1;
    memcpy(nodes[path - 2], reveal_path[AIMER_LOGN - depth], AIMER_SEED_SIZE);

    for (index = (1U << depth); index < (2U << depth); index += 4)
    {
      memcpy(&ctx, &ctx_, sizeof(hash_instance_x4));

      bufs[0][0] = (uint8_t)(rep_index);
      bufs[1][0] = (uint8_t)(rep_index);
      bufs[2][0] = (uint8_t)(rep_index);
      bufs[3][0] = (uint8_t)(rep_index);

      bufs[0][1] = (uint8_t)(index + 0);
      bufs[1][1] = (uint8_t)(index + 1);
      bufs[2][1] = (uint8_t)(index + 2);
      bufs[3][1] = (uint8_t)(index + 3);

      memcpy(bufs[0] + 2, nodes[index - 2], AIMER_SEED_SIZE);
      memcpy(bufs[1] + 2, nodes[index - 1], AIMER_SEED_SIZE);
      memcpy(bufs[2] + 2, nodes[index + 0], AIMER_SEED_SIZE);
      memcpy(bufs[3] + 2, nodes[index + 1], AIMER_SEED_SIZE);

      hash_update_x4(&ctx, in_ptrs, AIMER_SEED_SIZE + 2);
      hash_final_x4(&ctx);

      out_ptrs[0] = nodes[2 * index - 2];
      out_ptrs[1] = nodes[2 * index + 0];
      out_ptrs[2] = nodes[2 * index + 2];
      out_ptrs[3] = nodes[2 * index + 4];

      hash_squeeze_x4(&ctx, out_ptrs, AIMER_SEED_SIZE << 1);
    }
  }

  path = ((cover_index + AIMER_N) >> (AIMER_LOGN - depth)) ^ 1;
  memcpy(nodes[path - 2], reveal_path[AIMER_LOGN - depth], AIMER_SEED_SIZE);
}
