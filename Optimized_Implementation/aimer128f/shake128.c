#include <stdint.h>
#include "shake128.h"

const ALIGN2 uint64_t RC[200] = {
    0x0000000000000001,0x0000000000000001,0x0000000000000001,0x0000000000000001, 0x0000000000000001,0x0000000000000001,0x0000000000000001,0x0000000000000001,
    0x0000000000008082,0x0000000000008082,0x0000000000008082,0x0000000000008082, 0x0000000000008082,0x0000000000008082,0x0000000000008082,0x0000000000008082,
    0x800000000000808A,0x800000000000808A,0x800000000000808A,0x800000000000808A, 0x800000000000808A,0x800000000000808A,0x800000000000808A,0x800000000000808A,
    0x8000000080008000,0x8000000080008000,0x8000000080008000,0x8000000080008000, 0x8000000080008000,0x8000000080008000,0x8000000080008000,0x8000000080008000,
    0x000000000000808B,0x000000000000808B,0x000000000000808B,0x000000000000808B, 0x000000000000808B,0x000000000000808B,0x000000000000808B,0x000000000000808B,
    0x0000000080000001,0x0000000080000001,0x0000000080000001,0x0000000080000001, 0x0000000080000001,0x0000000080000001,0x0000000080000001,0x0000000080000001,
    0x8000000080008081,0x8000000080008081,0x8000000080008081,0x8000000080008081, 0x8000000080008081,0x8000000080008081,0x8000000080008081,0x8000000080008081,
    0x8000000000008009,0x8000000000008009,0x8000000000008009,0x8000000000008009, 0x8000000000008009,0x8000000000008009,0x8000000000008009,0x8000000000008009,
    0x000000000000008A,0x000000000000008A,0x000000000000008A,0x000000000000008A, 0x000000000000008A,0x000000000000008A,0x000000000000008A,0x000000000000008A,
    0x0000000000000088,0x0000000000000088,0x0000000000000088,0x0000000000000088, 0x0000000000000088,0x0000000000000088,0x0000000000000088,0x0000000000000088,
    0x0000000080008009,0x0000000080008009,0x0000000080008009,0x0000000080008009, 0x0000000080008009,0x0000000080008009,0x0000000080008009,0x0000000080008009,
    0x000000008000000A,0x000000008000000A,0x000000008000000A,0x000000008000000A, 0x000000008000000A,0x000000008000000A,0x000000008000000A,0x000000008000000A,
    0x000000008000808B,0x000000008000808B,0x000000008000808B,0x000000008000808B, 0x000000008000808B,0x000000008000808B,0x000000008000808B,0x000000008000808B,
    0x800000000000008B,0x800000000000008B,0x800000000000008B,0x800000000000008B, 0x800000000000008B,0x800000000000008B,0x800000000000008B,0x800000000000008B,
    0x8000000000008089,0x8000000000008089,0x8000000000008089,0x8000000000008089, 0x8000000000008089,0x8000000000008089,0x8000000000008089,0x8000000000008089,
    0x8000000000008003,0x8000000000008003,0x8000000000008003,0x8000000000008003, 0x8000000000008003,0x8000000000008003,0x8000000000008003,0x8000000000008003,
    0x8000000000008002,0x8000000000008002,0x8000000000008002,0x8000000000008002, 0x8000000000008002,0x8000000000008002,0x8000000000008002,0x8000000000008002,
    0x8000000000000080,0x8000000000000080,0x8000000000000080,0x8000000000000080, 0x8000000000000080,0x8000000000000080,0x8000000000000080,0x8000000000000080,
    0x000000000000800A,0x000000000000800A,0x000000000000800A,0x000000000000800A, 0x000000000000800A,0x000000000000800A,0x000000000000800A,0x000000000000800A,
    0x800000008000000A,0x800000008000000A,0x800000008000000A,0x800000008000000A, 0x800000008000000A,0x800000008000000A,0x800000008000000A,0x800000008000000A,
    0x8000000080008081,0x8000000080008081,0x8000000080008081,0x8000000080008081, 0x8000000080008081,0x8000000080008081,0x8000000080008081,0x8000000080008081,
    0x8000000000008080,0x8000000000008080,0x8000000000008080,0x8000000000008080, 0x8000000000008080,0x8000000000008080,0x8000000000008080,0x8000000000008080,
    0x0000000080000001,0x0000000080000001,0x0000000080000001,0x0000000080000001, 0x0000000080000001,0x0000000080000001,0x0000000080000001,0x0000000080000001,
    0x8000000080008008,0x8000000080008008,0x8000000080008008,0x8000000080008008, 0x8000000080008008,0x8000000080008008,0x8000000080008008,0x8000000080008008
  };

//   void print_512(__m512i x){
//     uint64_t c[8];
//     _mm512_store_si512((__m512i*)c,x);
//     printf("\n 0 = %8.16lx\t  1 = %8.16lx\t  2 = %8.16lx\t  3 = %8.16lx\t\n 4 = %8.16lx\t  5 = %8.16lx\t  6 = %8.16lx\t  7 = %8.16lx\t\n", c[0],c[1],c[2],c[3],c[4],c[5],c[6],c[7]);
// }

  /*organizar depois - função keccak*/

void keccakf(__m512i *s){
    __m512i aux[12];
    int i;
	for(i=0;i<24;i=i+2){
        aux[0] = _mm512_ternarylogic_epi64(s[0],s[5],s[10],0x96);
        aux[0] = _mm512_ternarylogic_epi64(aux[0],s[15],s[20],0x96);
        aux[2] = _mm512_ternarylogic_epi64(s[2],s[7],s[12],0x96);
        aux[2] = _mm512_ternarylogic_epi64(aux[2],s[17],s[22],0x96);
        aux[1] = _mm512_ternarylogic_epi64(s[1],s[6],s[11],0x96);
        aux[1] = _mm512_ternarylogic_epi64(aux[1],s[16],s[21],0x96);
        aux[3] = _mm512_ternarylogic_epi64(s[3],s[8],s[13],0x96);
        aux[3] = _mm512_ternarylogic_epi64(aux[3],s[18],s[23],0x96);
        aux[4] = _mm512_ternarylogic_epi64(s[4],s[9],s[14],0x96);
        aux[4] = _mm512_ternarylogic_epi64(aux[4],s[19],s[24],0x96);
        
        aux[6] = _mm512_rol_epi64(aux[2],0x1);/*aux[6] = d1*/
        aux[7] = _mm512_rol_epi64(aux[3],0x1);/*aux[7] = d2*/
        aux[9] = _mm512_rol_epi64(aux[0],0x1);/*aux[9] = d4*/
        aux[5] = _mm512_rol_epi64(aux[1],0x1);/*aux[5] = d0*/
        aux[8] = _mm512_rol_epi64(aux[4],0x1);/*aux[8] = d3*/
        
        s[0] = _mm512_ternarylogic_epi64(s[0],aux[4],aux[5],0x96);
        s[5] = _mm512_ternarylogic_epi64(s[5],aux[4],aux[5],0x96);
        s[10] = _mm512_ternarylogic_epi64(s[10],aux[4],aux[5],0x96);
        s[15] = _mm512_ternarylogic_epi64(s[15],aux[4],aux[5],0x96);
        s[20] = _mm512_ternarylogic_epi64(s[20],aux[4],aux[5],0x96);
        
        s[1] = _mm512_ternarylogic_epi64(s[1],aux[6],aux[0],0x96);
        s[6] = _mm512_ternarylogic_epi64(s[6],aux[6],aux[0],0x96);
        s[11] = _mm512_ternarylogic_epi64(s[11],aux[6],aux[0],0x96);
        s[16] = _mm512_ternarylogic_epi64(s[16],aux[6],aux[0],0x96);
        s[21] = _mm512_ternarylogic_epi64(s[21],aux[6],aux[0],0x96);
        
        s[2] = _mm512_ternarylogic_epi64(s[2],aux[7],aux[1],0x96);
        s[7] = _mm512_ternarylogic_epi64(s[7],aux[7],aux[1],0x96);
        s[12] = _mm512_ternarylogic_epi64(s[12],aux[7],aux[1],0x96);
        s[17] = _mm512_ternarylogic_epi64(s[17],aux[7],aux[1],0x96);
        s[22] = _mm512_ternarylogic_epi64(s[22],aux[7],aux[1],0x96);
        
        s[3] = _mm512_ternarylogic_epi64(s[3],aux[8],aux[2],0x96);
        s[8] = _mm512_ternarylogic_epi64(s[8],aux[8],aux[2],0x96);
        s[13] = _mm512_ternarylogic_epi64(s[13],aux[8],aux[2],0x96);
        s[18] = _mm512_ternarylogic_epi64(s[18],aux[8],aux[2],0x96);
        s[23] = _mm512_ternarylogic_epi64(s[23],aux[8],aux[2],0x96);
        
        s[4] = _mm512_ternarylogic_epi64(s[4],aux[9],aux[3],0x96);
        s[9] = _mm512_ternarylogic_epi64(s[9],aux[9],aux[3],0x96);
        s[14] = _mm512_ternarylogic_epi64(s[14],aux[9],aux[3],0x96);
        s[19] = _mm512_ternarylogic_epi64(s[19],aux[9],aux[3],0x96);
        s[24] = _mm512_ternarylogic_epi64(s[24],aux[9],aux[3],0x96);
        
        aux[1] = _mm512_rol_epi64(s[6],0x2C);
        aux[2] = _mm512_rol_epi64(s[12],0x2B);
        aux[3] = _mm512_rol_epi64(s[18],0x15);
        aux[4] = _mm512_rol_epi64(s[24],0x0E);
        
        
        s[18] = _mm512_ternarylogic_epi64(aux[4],s[0],aux[1],0xD2);/*s[18] = s4*/
        s[6] = _mm512_ternarylogic_epi64(aux[3],aux[4],s[0],0xD2);/*s[6] = s3*/
        s[0] = _mm512_ternarylogic_epi64(s[0],aux[1],aux[2],0xD2);/*s[0] = s0*/
        s[12] = _mm512_ternarylogic_epi64(aux[1],aux[2],aux[3],0xD2);/*s[12] = s1*/
        aux[6] = _mm512_loadu_si512((__m512i*)RC+i);
        s[24] = _mm512_ternarylogic_epi64(aux[2],aux[3],aux[4],0xD2);/*s[24] = s2*/
        
        s[0] = _mm512_xor_si512(s[0],aux[6]);
        
        
        aux[0] = _mm512_rol_epi64(s[3],0x1C);
        aux[1] = _mm512_rol_epi64(s[9],0x14);
        aux[2] = _mm512_rol_epi64(s[10],0x3);
        aux[3] = _mm512_rol_epi64(s[16],0x2D);
        aux[4] = _mm512_rol_epi64(s[22],0x3D);
        
        s[16] = _mm512_ternarylogic_epi64(aux[0],aux[1],aux[2],0xD2);/*s[16] = s5*/
        s[3]  = _mm512_ternarylogic_epi64(aux[1],aux[2],aux[3],0xD2);/*s[3] = s6*/
        s[10] = _mm512_ternarylogic_epi64(aux[2],aux[3],aux[4],0xD2);/*s[10] = s7*/
        s[22] = _mm512_ternarylogic_epi64(aux[3],aux[4],aux[0],0xD2);/*s[22] = s8*/
        s[9]  = _mm512_ternarylogic_epi64(aux[4],aux[0],aux[1],0xD2);/*s[9] = s9*/
        
        aux[0] = _mm512_rol_epi64(s[1],0x01);
        aux[1] = _mm512_rol_epi64(s[7],0x06);
        aux[2] = _mm512_rol_epi64(s[13],0x19);
        aux[3] = _mm512_rol_epi64(s[19],0x8);
        aux[4] = _mm512_rol_epi64(s[20],0x12);
        
        s[7]   = _mm512_ternarylogic_epi64(aux[0],aux[1],aux[2],0xD2);/*s[7] = s10*/
        s[19]  = _mm512_ternarylogic_epi64(aux[1],aux[2],aux[3],0xD2);/*s[19] = s11*/
        s[1]   = _mm512_ternarylogic_epi64(aux[2],aux[3],aux[4],0xD2);/*s[1] = s12*/
        s[13]  = _mm512_ternarylogic_epi64(aux[3],aux[4],aux[0],0xD2);/*s[13] = s13*/
        s[20]  = _mm512_ternarylogic_epi64(aux[4],aux[0],aux[1],0xD2);/*s[20] = s14*/
        
        aux[0] = _mm512_rol_epi64(s[4],0x1B);
        aux[1] = _mm512_rol_epi64(s[5],0x24);
        aux[2] = _mm512_rol_epi64(s[11],0xA);
        aux[3] = _mm512_rol_epi64(s[17],0xF);
        aux[4] = _mm512_rol_epi64(s[23],0x38);
        
        s[23] = _mm512_ternarylogic_epi64(aux[0],aux[1],aux[2],0xD2);/*s[23] = s15*/
        s[5]  = _mm512_ternarylogic_epi64(aux[1],aux[2],aux[3],0xD2);/*s[5] = s16*/
        s[17] = _mm512_ternarylogic_epi64(aux[2],aux[3],aux[4],0xD2);/*s[17] = s17*/
        s[4]  = _mm512_ternarylogic_epi64(aux[3],aux[4],aux[0],0xD2);/*s[4] = s18*/
        s[11] = _mm512_ternarylogic_epi64(aux[4],aux[0],aux[1],0xD2);/*s[11] = s19*/
        
        aux[0] = _mm512_rol_epi64(s[2],0x3E);
        aux[1] = _mm512_rol_epi64(s[8],0x37);
        aux[2] = _mm512_rol_epi64(s[14],0x27);
        aux[3] = _mm512_rol_epi64(s[15],0x29);
        aux[4] = _mm512_rol_epi64(s[21],0x2);
        
        s[14] = _mm512_ternarylogic_epi64(aux[0],aux[1],aux[2],0xD2);/*s[14] = s20*/
        s[21] = _mm512_ternarylogic_epi64(aux[1],aux[2],aux[3],0xD2);/*s[21] = s21*/
        s[8]  = _mm512_ternarylogic_epi64(aux[2],aux[3],aux[4],0xD2);/*s[8] = s22*/
        s[15] = _mm512_ternarylogic_epi64(aux[3],aux[4],aux[0],0xD2);/*s[15] = s23*/
        s[2]  = _mm512_ternarylogic_epi64(aux[4],aux[0],aux[1],0xD2);/*s[2] = s24*/
        
        /*new round*/
        
        aux[0] = _mm512_ternarylogic_epi64(s[0],s[16],s[7],0x96);
        aux[0] = _mm512_ternarylogic_epi64(aux[0],s[23],s[14],0x96);
        aux[2] = _mm512_ternarylogic_epi64(s[24],s[10],s[1],0x96);
        aux[2] = _mm512_ternarylogic_epi64(aux[2],s[17],s[8],0x96);
        aux[1] = _mm512_ternarylogic_epi64(s[12],s[3],s[19],0x96);
        aux[1] = _mm512_ternarylogic_epi64(aux[1],s[5],s[21],0x96);
        aux[3] = _mm512_ternarylogic_epi64(s[6],s[22],s[13],0x96);
        aux[3] = _mm512_ternarylogic_epi64(aux[3],s[4],s[15],0x96);
        aux[4] = _mm512_ternarylogic_epi64(s[18],s[9],s[20],0x96);
        aux[4] = _mm512_ternarylogic_epi64(aux[4],s[11],s[2],0x96);
        
        aux[6] = _mm512_rol_epi64(aux[2],0x01);/*aux[6] = d1*/
        aux[7] = _mm512_rol_epi64(aux[3],0x01);/*aux[7] = d2*/
        aux[9] = _mm512_rol_epi64(aux[0],0x01);/*aux[9] = d4*/
        aux[5] = _mm512_rol_epi64(aux[1],0x01);/*aux[5] = d0*/
        aux[8] = _mm512_rol_epi64(aux[4],0x01);/*aux[8] = d3*/
        
        
        s[0] = _mm512_ternarylogic_epi64(s[0],aux[5],aux[4],0x96);
        s[16] = _mm512_ternarylogic_epi64(s[16],aux[5],aux[4],0x96);
        s[7] = _mm512_ternarylogic_epi64(s[7],aux[5],aux[4],0x96);
        s[23] = _mm512_ternarylogic_epi64(s[23],aux[5],aux[4],0x96);
        s[14] = _mm512_ternarylogic_epi64(s[14],aux[5],aux[4],0x96);
        
        s[12] = _mm512_ternarylogic_epi64(s[12],aux[6],aux[0],0x96);
        s[3] = _mm512_ternarylogic_epi64(s[3],aux[6],aux[0],0x96);
        s[19] = _mm512_ternarylogic_epi64(s[19],aux[6],aux[0],0x96);
        s[5] = _mm512_ternarylogic_epi64(s[5],aux[6],aux[0],0x96);
        s[21] = _mm512_ternarylogic_epi64(s[21],aux[6],aux[0],0x96);
        
        s[24] = _mm512_ternarylogic_epi64(s[24],aux[7],aux[1],0x96);
        s[10] = _mm512_ternarylogic_epi64(s[10],aux[7],aux[1],0x96);
        s[1] = _mm512_ternarylogic_epi64(s[1],aux[7],aux[1],0x96);
        s[17] = _mm512_ternarylogic_epi64(s[17],aux[7],aux[1],0x96);
        s[8] = _mm512_ternarylogic_epi64(s[8],aux[7],aux[1],0x96);
        
        s[6] = _mm512_ternarylogic_epi64(s[6],aux[8],aux[2],0x96);
        s[22] = _mm512_ternarylogic_epi64(s[22],aux[8],aux[2],0x96);
        s[13] = _mm512_ternarylogic_epi64(s[13],aux[8],aux[2],0x96);
        s[4] = _mm512_ternarylogic_epi64(s[4],aux[8],aux[2],0x96);
        s[15] = _mm512_ternarylogic_epi64(s[15],aux[8],aux[2],0x96);
        
        s[18] = _mm512_ternarylogic_epi64(s[18],aux[9],aux[3],0x96);
        s[9] = _mm512_ternarylogic_epi64(s[9],aux[9],aux[3],0x96);
        s[20] = _mm512_ternarylogic_epi64(s[20],aux[9],aux[3],0x96);
        s[11] = _mm512_ternarylogic_epi64(s[11],aux[9],aux[3],0x96);
        s[2] = _mm512_ternarylogic_epi64(s[2],aux[9],aux[3],0x96);
        
        
        aux[1] = _mm512_rol_epi64(s[3],0x2C);
        aux[2] = _mm512_rol_epi64(s[1],0x2B);
        aux[3] = _mm512_rol_epi64(s[4],0x15);
        aux[4] = _mm512_rol_epi64(s[2],0x0E);
        
        s[3] = _mm512_ternarylogic_epi64(aux[3],aux[4],s[0],0xD2);/*s[3] = s3*/
        s[4] = _mm512_ternarylogic_epi64(aux[4],s[0],aux[1],0xD2);/*s[4] = s4*/
        s[0] = _mm512_ternarylogic_epi64(s[0],aux[1],aux[2],0xD2);/*s[0] = s0*/
        s[1] = _mm512_ternarylogic_epi64(aux[1],aux[2],aux[3],0xD2);/*s[1] = s1*/
        aux[6] = _mm512_loadu_si512((__m512i*)RC+i+1);
        s[2] = _mm512_ternarylogic_epi64(aux[2],aux[3],aux[4],0xD2);/*s[2] = s2*/
        
        s[0] = _mm512_xor_si512(s[0],aux[6]);
        
        aux[0] = _mm512_rol_epi64(s[6],0x1C);
        aux[1] = _mm512_rol_epi64(s[9],0x14);
        aux[2] = _mm512_rol_epi64(s[7],0x3);
        aux[3] = _mm512_rol_epi64(s[5],0x2D);
        aux[4] = _mm512_rol_epi64(s[8],0x3D);
        
        s[5] = _mm512_ternarylogic_epi64(aux[0],aux[1],aux[2],0xD2);/*s[5] = s5*/
        s[6] = _mm512_ternarylogic_epi64(aux[1],aux[2],aux[3],0xD2);/*s[6] = s6*/
        s[7] = _mm512_ternarylogic_epi64(aux[2],aux[3],aux[4],0xD2);/*s[7] = s7*/
        s[8] = _mm512_ternarylogic_epi64(aux[3],aux[4],aux[0],0xD2);/*s[8] = s8*/
        s[9] = _mm512_ternarylogic_epi64(aux[4],aux[0],aux[1],0xD2);/*s[9] = s9*/
        
        aux[0] = _mm512_rol_epi64(s[12],0x1);
        aux[1] = _mm512_rol_epi64(s[10],0x6);
        aux[2] = _mm512_rol_epi64(s[13],0x19);
        aux[3] = _mm512_rol_epi64(s[11],0x8);
        aux[4] = _mm512_rol_epi64(s[14],0x12);
        
        s[10]  = _mm512_ternarylogic_epi64(aux[0],aux[1],aux[2],0xD2);/*s[10] = s10*/
        s[11]  = _mm512_ternarylogic_epi64(aux[1],aux[2],aux[3],0xD2);/*s[11] = s11*/
        s[12]  = _mm512_ternarylogic_epi64(aux[2],aux[3],aux[4],0xD2);/*s[12] = s12*/
        s[13]  = _mm512_ternarylogic_epi64(aux[3],aux[4],aux[0],0xD2);/*s[13] = s13*/
        s[14]  = _mm512_ternarylogic_epi64(aux[4],aux[0],aux[1],0xD2);/*s[14] = s14*/
        
        aux[0] = _mm512_rol_epi64(s[18],0x1B);
        aux[1] = _mm512_rol_epi64(s[16],0x24);
        aux[2] = _mm512_rol_epi64(s[19],0xA);
        aux[3] = _mm512_rol_epi64(s[17],0xF);
        aux[4] = _mm512_rol_epi64(s[15],0x38);
        
        s[15] = _mm512_ternarylogic_epi64(aux[0],aux[1],aux[2],0xD2);/*s[15] = s15*/
        s[16]  = _mm512_ternarylogic_epi64(aux[1],aux[2],aux[3],0xD2);/*s[16] = s16*/
        s[17] = _mm512_ternarylogic_epi64(aux[2],aux[3],aux[4],0xD2);/*s[17] = s17*/
        s[18]  = _mm512_ternarylogic_epi64(aux[3],aux[4],aux[0],0xD2);/*s[18] = s18*/
        s[19] = _mm512_ternarylogic_epi64(aux[4],aux[0],aux[1],0xD2);/*s[19] = s19*/
        
        aux[0] = _mm512_rol_epi64(s[24],0x3E);
        aux[1] = _mm512_rol_epi64(s[22],0x37);
        aux[2] = _mm512_rol_epi64(s[20],0x27);
        aux[3] = _mm512_rol_epi64(s[23],0x29);
        aux[4] = _mm512_rol_epi64(s[21],0x2);
        
        s[20] = _mm512_ternarylogic_epi64(aux[0],aux[1],aux[2],0xD2);/*s[20] = s20*/
        s[21] = _mm512_ternarylogic_epi64(aux[1],aux[2],aux[3],0xD2);/*s[21] = s21*/
        s[22] = _mm512_ternarylogic_epi64(aux[2],aux[3],aux[4],0xD2);/*s[22] = s22*/
        s[23] = _mm512_ternarylogic_epi64(aux[3],aux[4],aux[0],0xD2);/*s[23] = s23*/
        s[24] = _mm512_ternarylogic_epi64(aux[4],aux[0],aux[1],0xD2);/*s[24] = s24*/
	}
}	

