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


  int num_blocks = len / BLOCK;
  if(len % BLOCK) num_blocks++;
  int num_states_align;

  int *all_states = (int *) _mm_malloc(sizeof(int)*num_states, 64);
  if(num_states % 16) num_states_align = 16 * (num_states / 16 +1);
  int *res_states = (int *) _mm_malloc(sizeof(int)*num_blocks*num_states_align, 64);
  for(int i=0;i<num_states;i++) all_states[i] = i;

  __m512i vsize = _mm512_set1_epi32(num_states);
  __m512i vspec = _mm512_set_epi32(0,0,0,0,0,0,0,0,0,6,5,4,3,2,1,0);


  struct timeval tv1, tv2;
  struct timezone tz1, tz2;
  gettimeofday(&tv1, &tz1);

  int submitted = 0;

#pragma omp parallel for num_threads(nth) 
  for(int j = 1; j < num_blocks; j++) {
    int start = j * BLOCK;
    int end = start + BLOCK < len ? start + BLOCK : len;


    __m512i vs = vspec;
    for(int i=start; i<end; i++) {
      __m512i vi = _mm512_set1_epi32(num_int[i]);
      vi = _mm512_add_epi32(vs, _mm512_mullo_epi32(vi, vsize));
      vs = _mm512_i32gather_epi32(vi, trans_table, 4);
    }
    _mm512_store_epi32(res_states+j*num_states_align, vs);
      //   for(;k<num_states;k++) {
      //    res_states[j*num_states_align+k] = trans_table[num_int[i]*num_states+res_states[j*num_states_align+k]];
      //  }
  }

  for(int i=0;i<BLOCK;i++) s = trans_table[num_int[i]*num_states+s];

  for(int j=1;j<num_blocks;j++) {
    s = res_states[j*num_states_align+s];
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
