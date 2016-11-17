#include <iostream>
#include <fstream>
#include <sstream>
#include <string>
#include <sys/time.h>
#include "immintrin.h"
using namespace std;

int num_states;
int *trans_table;

#define BLOCK 16384

void input_FSM(char *fsmfile) {
  ifstream fin(fsmfile);
  string line;
  getline(fin, line);
  stringstream sin1(line);
  sin1 >> num_states;

  trans_table = (int *)_mm_malloc(sizeof(int)*num_states*2, 64);

  int s, c, e;

  while(getline(fin, line)) {
    stringstream sin(line);
    sin >> s >> c >> e;
    trans_table[c*num_states+s] = e;
  }
}


bool check_div(string num, int nth) {
  int *num_int = (int *)_mm_malloc(sizeof(int)*num.length(), 64);
  for(int i=0;i<num.length();i++) {
    num_int[i] = num[i] - '0';
  }

  int len = num.length();

  int s = 0;

 int *spec_temp = (int *)_mm_malloc(sizeof(int)*16, 64);

 for(int i=0;i<16;i++) spec_temp[i] = 0;
 for(int i=0;i<8 && i<num_states;i++) {
  spec_temp[i] = i;
 }
 for(int i=8;i<16 && i-8 < num_states;i++) {
  spec_temp[i] = i-8;
 }

  __m512i vspec = _mm512_load_epi32(spec_temp); 

 for(int i=0;i<16;i++) spec_temp[i] = 0;
 for(int i=0;i<16 && i<num_states;i++) {
  spec_temp[i] = i;
 }
  __m512i vspec1 = _mm512_load_epi32(spec_temp); 

  int num_blocks = len / BLOCK;
  __m512i *res_states = (__m512i *) _mm_malloc(sizeof(__m512i)*num_blocks, 64);
  __m512i vperm = _mm512_set_epi32(7,6,5,4,3,2,1,0,15,14,13,12,11,10,9,8);

  int submitted = 0;


  struct timezone tz1, tz2;
  struct timeval tv1, tv2;
  gettimeofday(&tv1, &tz1);
  int nb = num_blocks % 2 == 0 ? num_blocks-1:num_blocks;
#pragma omp parallel for num_threads(nth) 
  for(int j = 1; j < nb; j+=2) {
    int start = j * BLOCK;
    int end = start + BLOCK;
    __m512i vs = vspec;

    for(int i=start; i<end; i+=8) {
      __m512i a = _mm512_set1_epi32(num_int[i]);
      __m512i b = _mm512_set1_epi32(num_int[i+BLOCK]);
      __m512i vi = _mm512_mask_permutevar_epi32(a, 0xFF00, vperm, b);
      vi = _mm512_mullo_epi32(vi, _mm512_set1_epi32(num_states));
      vi = _mm512_add_epi32(vs, vi);
      vs = _mm512_i32gather_epi32(vi, trans_table, 4);

      __m512i a1 = _mm512_set1_epi32(num_int[i+1]);
      __m512i b1 = _mm512_set1_epi32(num_int[i+1+BLOCK]);
      __m512i vi1 = _mm512_mask_permutevar_epi32(a1, 0xFF00, vperm, b1);
      vi1 = _mm512_mullo_epi32(vi1, _mm512_set1_epi32(num_states));
      vi1 = _mm512_add_epi32(vs, vi1);
      vs = _mm512_i32gather_epi32(vi1, trans_table, 4);

      __m512i a2 = _mm512_set1_epi32(num_int[i+2]);
      __m512i b2 = _mm512_set1_epi32(num_int[i+2+BLOCK]);
      __m512i vi2 = _mm512_mask_permutevar_epi32(a2, 0xFF00, vperm, b2);
      vi2 = _mm512_mullo_epi32(vi2, _mm512_set1_epi32(num_states));
      vi2 = _mm512_add_epi32(vs, vi2);
      vs = _mm512_i32gather_epi32(vi2, trans_table, 4);

      __m512i a3 = _mm512_set1_epi32(num_int[i+3]);
      __m512i b3 = _mm512_set1_epi32(num_int[i+3+BLOCK]);
      __m512i vi3 = _mm512_mask_permutevar_epi32(a3, 0xFF00, vperm, b3);
      vi3 = _mm512_mullo_epi32(vi3, _mm512_set1_epi32(num_states));
      vi3 = _mm512_add_epi32(vs, vi3);
      vs = _mm512_i32gather_epi32(vi3, trans_table, 4);

      __m512i a4 = _mm512_set1_epi32(num_int[i+4]);
      __m512i b4 = _mm512_set1_epi32(num_int[i+4+BLOCK]);
      __m512i vi4 = _mm512_mask_permutevar_epi32(a4, 0xFF00, vperm, b4);
      vi4 = _mm512_mullo_epi32(vi4, _mm512_set1_epi32(num_states));
      vi4 = _mm512_add_epi32(vs, vi4);
      vs = _mm512_i32gather_epi32(vi4, trans_table, 4);

      __m512i a5 = _mm512_set1_epi32(num_int[i+5]);
      __m512i b5 = _mm512_set1_epi32(num_int[i+5+BLOCK]);
      __m512i vi5 = _mm512_mask_permutevar_epi32(a5, 0xFF00, vperm, b5);
      vi5 = _mm512_mullo_epi32(vi5, _mm512_set1_epi32(num_states));
      vi5 = _mm512_add_epi32(vs, vi5);
      vs = _mm512_i32gather_epi32(vi5, trans_table, 4);

      __m512i a6 = _mm512_set1_epi32(num_int[i+6]);
      __m512i b6 = _mm512_set1_epi32(num_int[i+6+BLOCK]);
      __m512i vi6 = _mm512_mask_permutevar_epi32(a6, 0xFF00, vperm, b6);
      vi6 = _mm512_mullo_epi32(vi6, _mm512_set1_epi32(num_states));
      vi6 = _mm512_add_epi32(vs, vi6);
      vs = _mm512_i32gather_epi32(vi6, trans_table, 4);

      __m512i a7 = _mm512_set1_epi32(num_int[i+7]);
      __m512i b7 = _mm512_set1_epi32(num_int[i+7+BLOCK]);
      __m512i vi7 = _mm512_mask_permutevar_epi32(a7, 0xFF00, vperm, b7);
      vi7 = _mm512_mullo_epi32(vi7, _mm512_set1_epi32(num_states));
      vi7 = _mm512_add_epi32(vs, vi7);
      vs = _mm512_i32gather_epi32(vi7, trans_table, 4);


    }

    res_states[j] = vs;
    res_states[j+1] = _mm512_permutevar_epi32(vperm, vs);
  }

  if(num_blocks % 2 == 0) {
    int start = (num_blocks-1) * BLOCK;
    int end = start + BLOCK;
    __m512i vs = vspec1;
    for(int i=start; i<end; i++) {
      __m512i vi = _mm512_set1_epi32(num_int[i]);
      vi = _mm512_mullo_epi32(vi, _mm512_set1_epi32(num_states));
      vi = _mm512_add_epi32(vs, vi);
      vs = _mm512_i32gather_epi32(vi, trans_table, 4);
    }
    res_states[num_blocks-1] = vs;
  }

  int cur = 0;
  int ppp = 0;


  while(cur < len) {
    int end = cur + BLOCK < len ? cur + BLOCK : len;
    for(int i=cur;i<end;i++) {
      s = trans_table[s+(num_int[i])*num_states];
    }
    ppp++;
    submitted++;
    cur += BLOCK;
    while(submitted < num_blocks) {
      if(s < 8) {
        s = _mm512_mask_reduce_add_epi32(1<<(s+1), res_states[submitted]); 
        submitted++;
        cur += BLOCK;
      }
      else {
        break; 
      }
    }
  }
  gettimeofday(&tv2, &tz2);
  cout << "used time: " << tv2.tv_usec - tv1.tv_usec + 1000000 * (tv2.tv_sec - tv1.tv_sec) << endl;
  if(s==0)return true;
  return false;
}

int main(int argc, char *argv[])
{
  input_FSM("datasets/div.fsm");

  ifstream fin(argv[1]);
  string num;
  getline(fin, num);

  int nth = atoi(argv[2]);
  bool divisable = check_div(num, nth);
  return 0;
}
