#include <iostream>
#include <fstream>
#include <sstream>
#include <string>
#include <sys/time.h>
#include "immintrin.h"
using namespace std;

#define BLOCK 8
__mmask16 least_one[65536];
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
    for(int j=1;j<16;j++)tmp_states[j] = (j-1) % 7;
    spec_states[i] = _mm512_load_epi32(tmp_states);
  }

  vind = (__m512i *)_mm_malloc(sizeof(__m512i)*BLOCK, 64);
  for(int i=0;i<BLOCK;i++) {
    tmp_states[0] = i;
    for(int j=1;j<16;j++)tmp_states[j] = i+BLOCK;
    vind[i] = _mm512_load_epi32(tmp_states);
  }

  for(int i=2;i<65536;i++) {
    int c = 0;
    int t = i & 0xFFFE;
    while(t % 2 == 0) {
      t /= 2;
      c++;
    }
    least_one[i] = (1 << c);
  }


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

  __m512i vs;
  __m512i vss;
  int *vss_buffer = (int *)_mm_malloc(sizeof(int)*16, 64);


  struct timeval tv1, tv2;
  struct timezone tz1, tz2;
  gettimeofday(&tv1, &tz1);

  while(cur + 2*BLOCK < num.length()) {
    vs = spec_states[s];
    vss = vs;

    int i= cur;
    __m512i vp1 = _mm512_i32gather_epi32(vind[i-cur], num_int+cur, 4);
    vp1 = _mm512_mullo_epi32(vp1, _mm512_set1_epi32(num_states));
    vp1 = _mm512_add_epi32(vss, vp1);
    vss = _mm512_i32gather_epi32(vp1, trans_table, 4);

    __m512i vp2 = _mm512_i32gather_epi32(vind[i+1-cur], num_int+cur, 4);
    vp2 = _mm512_mullo_epi32(vp2, _mm512_set1_epi32(num_states));
    vp2 = _mm512_add_epi32(vss, vp2);
    vss = _mm512_i32gather_epi32(vp2, trans_table, 4);

    __m512i vp3 = _mm512_i32gather_epi32(vind[i+2-cur], num_int+cur, 4);
    vp3 = _mm512_mullo_epi32(vp3, _mm512_set1_epi32(num_states));
    vp3 = _mm512_add_epi32(vss, vp3);
    vss = _mm512_i32gather_epi32(vp3, trans_table, 4);

    __m512i vp4 = _mm512_i32gather_epi32(vind[i+3-cur], num_int+cur, 4);
    vp4 = _mm512_mullo_epi32(vp4, _mm512_set1_epi32(num_states));
    vp4 = _mm512_add_epi32(vss, vp4);
    vss = _mm512_i32gather_epi32(vp4, trans_table, 4);

    __m512i vp5 = _mm512_i32gather_epi32(vind[i+4-cur], num_int+cur, 4);
    vp5 = _mm512_mullo_epi32(vp5, _mm512_set1_epi32(num_states));
    vp5 = _mm512_add_epi32(vss, vp5);
    vss = _mm512_i32gather_epi32(vp5, trans_table, 4);

    __m512i vp6 = _mm512_i32gather_epi32(vind[i+5-cur], num_int+cur, 4);
    vp6 = _mm512_mullo_epi32(vp6, _mm512_set1_epi32(num_states));
    vp6 = _mm512_add_epi32(vss, vp6);
    vss = _mm512_i32gather_epi32(vp6, trans_table, 4);

    __m512i vp7 = _mm512_i32gather_epi32(vind[i+6-cur], num_int+cur, 4);
    vp7 = _mm512_mullo_epi32(vp7, _mm512_set1_epi32(num_states));
    vp7 = _mm512_add_epi32(vss, vp7);
    vss = _mm512_i32gather_epi32(vp7, trans_table, 4);

    __m512i vp8 = _mm512_i32gather_epi32(vind[i+7-cur], num_int+cur, 4);
    vp8 = _mm512_mullo_epi32(vp8, _mm512_set1_epi32(num_states));
    vp8 = _mm512_add_epi32(vss, vp8);
    vss = _mm512_i32gather_epi32(vp8, trans_table, 4);


    _mm512_store_epi32(vss_buffer, vss);
    ppp+=2;
    qqq++;
    cur += BLOCK;
    if(0 == vss_buffer[0]) {
      s = vss_buffer[2];
      cur += BLOCK;
      qqq++;
    }
  }

  for(int i=cur;i<num.length();i++) {
    s = trans_table[num_int[i]*num_states+s];
  }
  gettimeofday(&tv2, &tz2);
  cout << "used time: " << tv2.tv_usec - tv1.tv_usec + 1000000 * (tv2.tv_sec - tv1.tv_sec) << endl;
  cout << "spec success rate: " << qqq / ppp << endl;
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
  cout << "divisable: " << divisable << endl;
  return 0;
}
