#include <iostream>
#include <fstream>
#include <sstream>
#include <map>
#include <sys/time.h>
#include <stdlib.h>
#include "immintrin.h"

using namespace std;

int num_states;
int num_alphs;
int *trans_table;
map<char, int> lookup_table;
int cur_in = 0;

#define BLOCK 4096

void input_rex(char *filename)
{
 ifstream fin(filename); 
 string line;
 getline(fin, line);
 stringstream sin(line);
 sin >> num_states >> num_alphs;
 trans_table = (int *)_mm_malloc(sizeof(int)*num_states*(num_alphs+1), 64);
 for(int i=0;i<=num_alphs;i++) {
  for(int j=0;j<num_states;j++) {
    trans_table[i*num_states+j] = 0;
  }
 }

 int ss;
 getline(fin, line);
 stringstream sin1(line);
// while(sin1 >> ss) {
//   trans_table[num_alphs*num_states+ss] = 1;
// }

 while(getline(fin, line)) {
  stringstream sin2(line);
  int s, e;
  char in;
  sin2 >> s >> in >> e;

  if(in=='^') {
    for(int i=0;i<=num_alphs;i++) trans_table[i*num_states+s] = e;
  } else {
    map<char, int>::iterator it = lookup_table.find(in);
    if(it != lookup_table.end())
      trans_table[(it->second)*num_states+s] = e;
    else {
      lookup_table[in] = cur_in;
      trans_table[cur_in*num_states+s] = e;
      cur_in++;
    }
  }
 }
}

void print_rex() 
{
  for(int i=0;i<=num_alphs;i++) {
    for(int j=0;j<num_states;j++) {
      cout << trans_table[i*num_states+j]<< " ";
    }
    cout << endl;
  }
}

void match(char *filename, int nth) 
{
  ifstream fin(filename);
  char c;
  int s = 0;
  int count = 0;
  string text;
  string line;

  while(getline(fin, line)) {
    text += line;
  }

  for(int i=0;i<text.length();i++) {
    if(text[i]==' ')text[i] = '|';
  }

  int *text_int = (int *)_mm_malloc(sizeof(int)*text.length(), 64);

  for(int i=0;i<text.length();i++) {
    map<char, int>::iterator it = lookup_table.find(text[i]);
    if(it!=lookup_table.end()) {
      text_int[i] = it->second;
    } else {
      text_int[i] = cur_in;
    }
  }


 int *spec_temp = (int *)_mm_malloc(sizeof(int)*16, 64);
 int i = 0;
 for(;i<16 && i<num_states;i++) {
  spec_temp[i] = i;
 }
 for(;i<16;i++) {
  spec_temp[i] = 0;
 }

  __m512i vspec = _mm512_load_epi32(spec_temp); 

  int num_blocks = text.length() / BLOCK;
  if(text.length() %BLOCK) num_blocks++;
  int num_states_align;

  int *all_states = (int *) _mm_malloc(sizeof(int)*num_states, 64);

  for(int i=0;i<num_states;i++) all_states[i] = i;
  if(num_states % 16) num_states_align = 16 * (num_states / 16 +1);
  int *res_states = (int *) _mm_malloc(sizeof(int)*num_blocks*num_states_align, 64);

  __m512i vsize = _mm512_set1_epi32(num_states);
  struct timezone tz1, tz2;
  struct timeval tv1, tv2;
  gettimeofday(&tv1, &tz1);

  int submitted = 0;

#pragma omp parallel for num_threads(nth) 
  for(int j = 1; j < num_blocks; j++) {
    int start = j * BLOCK;
    int end = start + BLOCK < text.length() ? start + BLOCK : text.length();

    for(int i=0;i<num_states;i++) res_states[j*num_states_align+i] = all_states[i];

    for(int i=start; i<start+32; i++) {
      int k=0;
      for(;k<num_states/16*16;k+=16) {
        __m512i vs = _mm512_load_epi32(res_states+j*num_states_align+k);
        __m512i vi = _mm512_set1_epi32(text_int[i]);
        vi = _mm512_add_epi32(vs, _mm512_mullo_epi32(vi, vsize));
        vs = _mm512_i32gather_epi32(vi, trans_table, 4);
        _mm512_store_epi32(res_states+j*num_states_align+k, vs);
      }
      for(;k<num_states;k++) {
        res_states[j*num_states_align+k] = trans_table[text_int[i]*num_states+res_states[j*num_states_align+k]];
      }
    }

    int size_of_states1 = 1;

    for(int i=32+start;i<end;i++) {
      int k=0;
      for(;k<size_of_states1;k++) {
        res_states[j*num_states_align+k] = trans_table[text_int[i]*num_states+res_states[j*num_states_align+k]];
      }
    }


  }


  int cur = 0;
  int ppp = 0;

  for(int i=0;i<BLOCK;i++) s = trans_table[text_int[i]*num_states+s];

  for(int j=1;j<num_blocks;j++) {
    s = res_states[j*num_states_align+0];
  }



  gettimeofday(&tv2, &tz2);
  cout << "final state: " << s << endl;
  cout << "time used: " << tv2.tv_usec - tv1.tv_usec + 1000000 * (tv2.tv_sec - tv1.tv_sec) << "ms"<<endl;

}


int main(int argc, char *argv[])
{
  int nth = atoi(argv[3]);
  input_rex(argv[1]);
  match(argv[2], nth);


  return 0;
}
