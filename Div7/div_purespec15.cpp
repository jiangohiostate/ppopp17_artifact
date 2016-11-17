#include <iostream>
#include <fstream>
#include <sstream>
#include <string>
#include <sys/time.h>
#include "immintrin.h"
using namespace std;

#define BLOCK 8
__mmask16 least_one[65536];
int least_one_pos[65536];
__m512i *spec_states;
__m512i *vind;

int num_states;
int *trans_table;

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
  

  spec_states = (__m512i *)_mm_malloc(sizeof(__m512i)*num_states, 64);
  int *tmp_states = (int *)_mm_malloc(sizeof(int)*16, 64);
  for(int i=0;i<num_states;i++) {
    tmp_states[0] = i;
    for(int j=1;j<16;j++)tmp_states[j] = (j-1) % 5;
    spec_states[i] = _mm512_load_epi32(tmp_states);
  }

  vind = (__m512i *)_mm_malloc(sizeof(__m512i)*BLOCK, 64);
  for(int i=0;i<BLOCK;i++) {
    tmp_states[0] = i;
    for(int j=1;j<16;j++)tmp_states[j] = ((j-1) / 5 + 1) * BLOCK;
    vind[i] = _mm512_load_epi32(tmp_states);
  }
  for(int i=1;i<65536;i++) {
    int c = 0;
    int t = i;
    while(t % 2 == 0) {
      t /= 2;
      c++;
    }
    least_one[i] = (1 << c);
    least_one_pos[i] = c+1;
  }

  least_one_pos[0] = 16;
}


bool check_div(string num) {
  int *num_int = (int *)_mm_malloc(sizeof(int)*num.length(), 64);
  for(int i=0;i<num.length();i++) {
    num_int[i] = num[i] - '0';
  }

  int s = 0;

  int cur = 0;

  float ppp = 0;
  float qqq = 0;

  __m512i vp;
  __m512i vs;
  __m512i vss;
  int *vss_buffer = (int *)_mm_malloc(sizeof(int)*16, 64);

  __m512i vzero = _mm512_set1_epi32(0);

  struct timeval tv1, tv2;
  struct timezone tz1, tz2;
  gettimeofday(&tv1, &tz1);

  while(cur + 4*BLOCK < num.length()) {
    vs = _mm512_set1_epi32(s);
    vs = _mm512_mask_mov_epi32(vs, 0xFFFE, vzero);
    vss = vs;

    for(int i=cur;i<cur+BLOCK;i++) {
      vp = _mm512_i32gather_epi32(vind[i-cur], num_int+cur, 4);
      vp = _mm512_mullo_epi32(vp, _mm512_set1_epi32(num_states));
      vp = _mm512_add_epi32(vss, vp);
      vss = _mm512_i32gather_epi32(vp, trans_table, 4);
    }
    ppp+=16;
    qqq++;
    __mmask16 m = _mm512_cmpneq_epi32_mask(vzero, vss);
    if(m==0) {
      s = 0;
      cur+=16*BLOCK;
      qqq+=15;
    } else {
      s = _mm512_mask_reduce_add_epi32(least_one[m], vss);
      cur += least_one_pos[m] * BLOCK;
      qqq+=least_one_pos[m];
    }
  }

  for(int i=cur;i<num.length();i++) {
    s = trans_table[num_int[i]*num_states+s];
  }
  gettimeofday(&tv2, &tz2);
  cout << "used time: " << tv2.tv_usec - tv1.tv_usec + 1000000 * (tv2.tv_sec - tv1.tv_sec) << endl;
  cout << "spec succ rate: " << qqq/ppp << endl;
  if(s==0)return true;
  return false;
}

int main(int argc, char *argv[])
{
  input_FSM("datasets/div.fsm");
  
  ifstream fin(argv[1]);
  string num;
  getline(fin, num);

  bool divisable = check_div(num);
  return 0;
}
