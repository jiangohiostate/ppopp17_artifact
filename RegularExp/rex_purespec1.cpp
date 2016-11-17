#include <iostream>
#include <fstream>
#include <sstream>
#include <map>
#include <sys/time.h>
#include <limits.h>
#include "immintrin.h"
using namespace std;

#define BLOCK 8

int num_states;
int num_alphs;
int *trans_table;
int *spec_table_temp;
__m512i *vind;
map<char, int> lookup_table;
int cur_in = 0;
__m512i *spec_mm_table;


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

 spec_mm_table = (__m512i *)_mm_malloc(sizeof(__m512i)*num_states, 64);

 int *spec_temp = (int *)_mm_malloc(sizeof(int)*16, 64);
 spec_table_temp = (int *)_mm_malloc(sizeof(int)*(num_alphs+1)*(num_alphs+1), 64);

 int *all_states = (int *)_mm_malloc(sizeof(int)*num_states, 64);
 int *to_states = (int *)_mm_malloc(sizeof(int)*num_states, 64);
 int *count = new int[num_states];
 int max, p;


 for(int i1=0;i1<=num_alphs;i1++) {
   for(int i2=0;i2<=num_alphs;i2++) {
     for(int k=0;k<num_states;k++) {
       to_states[k] = trans_table[trans_table[i1*num_states+k]+i2*num_states]; 
     }
     for(int k=0;k<num_states;k++) {
       count[k] = 0; 
     }

     for(int k=0;k<num_states;k++) {
       count[to_states[k]]++; 
     }

     max = INT_MIN;
     for(int k=0;k<num_states;k++) {
       if(count[k] > max) {
         max = count[k];
         p = k;
       }
     }
     spec_table_temp[i1*(num_alphs+1)+i2] = p;
   }
 }




 for(int i=0;i<num_states;i++) {
   spec_temp[0] = i;
   int j = 1;
   for(;j<16 && j-1<num_states;j++) {
     spec_temp[j] = j-1;
   }
   for(;j<16;j++) {
     spec_temp[j] = 0;
   }
   spec_mm_table[i] = _mm512_load_epi32(spec_temp);
 }


 int *tmp_states = (int *)_mm_malloc(sizeof(int)*16, 64);
 vind = (__m512i *)_mm_malloc(sizeof(__m512i)*BLOCK, 64);
 for(int i=0;i<BLOCK;i++) {
   tmp_states[0] = i;
   for(int j=1;j<16;j++)tmp_states[j] = i+BLOCK;
   vind[i] = _mm512_load_epi32(tmp_states);
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

void match(char *filename) 
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

  int *text_int = (int *)malloc(sizeof(int) * text.length());
  for(int i=0;i<text.length();i++) {
    map<char, int>::iterator it = lookup_table.find(text[i]);
    if(it!=lookup_table.end()) {
      text_int[i] = it->second;
    } else {
      text_int[i] = num_alphs;
    }
  }

  __m512i vsize = _mm512_set1_epi32(num_states);
  int *tmp_states = (int *) _mm_malloc(sizeof(int)*16, 64);
  __m512i vperm = _mm512_set_epi32(7,6,5,4,3,2,1,0,15,14,13,12,11,10,9,8);
  __m512i vperm1 = _mm512_set_epi32(14,13,12,11,10,9,8,7,6,5,4,3,2,1,0,15);
  for(int i=0;i<16;i++) tmp_states[i] = 0;
  for(int i=0;i<16&&i<num_states;i++) tmp_states[i] = i;
  __m512i b = _mm512_load_epi32(tmp_states);

  float ppp=0,qqq=0;
  int cur = 0;
  __m512i vs;
  __m512i vss;
  int *vss_buffer = (int *)_mm_malloc(sizeof(int)*16, 64);
  struct timeval tv1, tv2;
  struct timezone tz1, tz2;
  gettimeofday(&tv1, &tz1);

  while(cur + 2*BLOCK < text.length()) {
    vs = _mm512_set1_epi32(s);
    int ss = spec_table_temp[text_int[cur+BLOCK-2]*(num_alphs+1)+text_int[cur+BLOCK-1]];
    vs = _mm512_mask_mov_epi32(vs, 0x0002, _mm512_set1_epi32(ss));
    vss = vs;

    for(int i=cur;i<cur+BLOCK;i++) {
      __m512i vp1 = _mm512_i32gather_epi32(vind[i-cur], text_int+cur, 4);
      vp1 = _mm512_mullo_epi32(vp1, vsize);
      vp1 = _mm512_add_epi32(vss, vp1);
      vss = _mm512_i32gather_epi32(vp1, trans_table, 4);
    }

    _mm512_store_epi32(vss_buffer, vss);
    s = vss_buffer[0];
    cur += BLOCK;
    ppp+=2;
qqq++;
    if(s==ss) {
      qqq++;
      s = vss_buffer[1];
      cur += BLOCK;
    }
  }

  for(int i=cur;i<text.length();i++) {
      s = trans_table[text_int[i]*num_states+s];
  }

  gettimeofday(&tv2, &tz2);
  cout << "time used: " << tv2.tv_usec - tv1.tv_usec + 1000000*(tv2.tv_sec - tv1.tv_sec) << "ms"<< endl;
  cout << "final state: " << s << endl;
  cout << "success rate: " << qqq / ppp << endl;
}


int main(int argc, char *argv[])
{
  input_rex(argv[1]);
  match(argv[2]);
  return 0;
}
