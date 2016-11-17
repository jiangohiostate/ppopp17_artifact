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
map<char, int> lookup_table;
int cur_in = 0;
__m512i *spec_mm_table;
__m512i *input_mm_table;

__mmask16 least_one[65536];
int least_one_pos[65536];
__m512i *vind;

void print_vector(__m512i v) {
  for(int i=0;i<16;i++) {
    cout << *((int *)&v + i) << " ";
  }
  cout << endl;
}

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
 for(int i=0;i<num_states;i++) {
   spec_temp[0] = i;
   int j = 1;
   for(;j<4 && j-1<num_states;j++) {
     spec_temp[j] = j-1;
   }
   for(;j<4;j++) {
     spec_temp[j] = 0;
   }
   spec_temp[4] = 0;
   spec_temp[5] = 1;
   spec_temp[6] = 0;
   spec_temp[7] = 1;
   spec_temp[8] = 0;
   spec_temp[9] = 1;
   spec_temp[10] = 0;
   spec_temp[11] = 1;
   spec_temp[12] = 0;
   spec_temp[13] = 1;
   spec_temp[14] = 0;
   spec_temp[15] = 1;

   spec_mm_table[i] = _mm512_load_epi32(spec_temp);
 }

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






  input_mm_table = (__m512i *)_mm_malloc(sizeof(__m512i)*num_alphs*num_alphs, 64);
 int *input_temp = (int *)_mm_malloc(sizeof(int)*16, 64);
  for(int i=0;i<num_alphs;i++) {
    for(int j=0;j<num_alphs;j++) {
      input_temp[0] = i*num_states;
      for(int k=1;k<16;k++) {
        input_temp[k] = j*num_states;
      }
      input_mm_table[i*num_alphs+j] = _mm512_load_epi32(input_temp);
    }
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
  vind = (__m512i *)_mm_malloc(sizeof(__m512i)*BLOCK, 64);
  int *tmp_states = (int *)_mm_malloc(sizeof(int)*16, 64);
  for(int i=0;i<BLOCK;i++) {
    tmp_states[0] = i;
    for(int j=1;j<4;j++)tmp_states[j] = i+BLOCK;
    for(int j=4;j<16;j++)tmp_states[j] = i+j/2*BLOCK;
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
  int *text_int = (int *)_mm_malloc(sizeof(int)*text.length(), 64);

  for(int i=0;i<text.length();i++) {
    map<char, int>::iterator it = lookup_table.find(text[i]);
    if(it!=lookup_table.end()) {
      text_int[i] = it->second;
    } else {
      text_int[i] = cur_in;
    }
  }

  float ppp = 0, qqq = 0;
  int cur = 0;
  __m512i vp;
  __m512i vs;
  __m512i vss;
  __m512i vzero = _mm512_set1_epi32(0);
  __m512i vtwo = _mm512_set1_epi32(2);
  __m512i vsize = _mm512_set1_epi32(num_states);
  int *input_temp = (int *)_mm_malloc(sizeof(int)*16, 64);
  int *vss_buffer = (int *)_mm_malloc(sizeof(int)*16, 64);

  __m512i vblock = _mm512_set_epi32(15*BLOCK, 14*BLOCK, 13*BLOCK, 12*BLOCK, 11*BLOCK, 10*BLOCK, 9*BLOCK, 8*BLOCK, 7*BLOCK, 6*BLOCK, 5*BLOCK, 4*BLOCK, 3*BLOCK, 2*BLOCK, BLOCK, 0);
  struct timeval tv1, tv2;

  struct timezone tz1, tz2;
  gettimeofday(&tv1, &tz1);
  
  while(cur + 16*BLOCK < text.length()) {
    vs = _mm512_set1_epi32(s);
   __m512i vss1 = _mm512_set1_epi32(spec_table_temp[text_int[cur+BLOCK-2]*(num_alphs+1)+text_int[cur+BLOCK-1]]);
    vs = _mm512_mask_mov_epi32(vs, 0xFFFE, vss1);
    vss = vs;

    for(int i=cur;i<cur+BLOCK;i+=8) {
      __m512i vp1 = _mm512_i32gather_epi32(_mm512_add_epi32(vblock, _mm512_set1_epi32(i)), &(text_int[0]), 4);
      vp1 = _mm512_mullo_epi32(vp1, vsize);
      vp1 = _mm512_add_epi32(vss, vp1);
      vss = _mm512_i32gather_epi32(vp1, trans_table, 4);

     __m512i vp2 = _mm512_i32gather_epi32(_mm512_add_epi32(vblock, _mm512_set1_epi32(i+1)), &(text_int[0]), 4);
      vp2 = _mm512_mullo_epi32(vp2, vsize);
      vp2 = _mm512_add_epi32(vss, vp2);
      vss = _mm512_i32gather_epi32(vp2, trans_table, 4);

      __m512i vp3 = _mm512_i32gather_epi32(_mm512_add_epi32(vblock, _mm512_set1_epi32(i+2)), &(text_int[0]), 4);
      vp3 = _mm512_mullo_epi32(vp3, vsize);
      vp3 = _mm512_add_epi32(vss, vp3);
      vss = _mm512_i32gather_epi32(vp3, trans_table, 4);

      __m512i vp4 = _mm512_i32gather_epi32(_mm512_add_epi32(vblock, _mm512_set1_epi32(i+3)), &(text_int[0]), 4);
      vp4 = _mm512_mullo_epi32(vp4, vsize);
      vp4 = _mm512_add_epi32(vss, vp4);
      vss = _mm512_i32gather_epi32(vp4, trans_table, 4);

      __m512i vp5 = _mm512_i32gather_epi32(_mm512_add_epi32(vblock, _mm512_set1_epi32(i+4)), &(text_int[0]), 4);
      vp5 = _mm512_mullo_epi32(vp5, vsize);
      vp5 = _mm512_add_epi32(vss, vp5);
      vss = _mm512_i32gather_epi32(vp5, trans_table, 4);

      __m512i vp6 = _mm512_i32gather_epi32(_mm512_add_epi32(vblock, _mm512_set1_epi32(i+5)), &(text_int[0]), 4);
      vp6 = _mm512_mullo_epi32(vp6, vsize);
      vp6 = _mm512_add_epi32(vss, vp6);
      vss = _mm512_i32gather_epi32(vp6, trans_table, 4);

      __m512i vp7 = _mm512_i32gather_epi32(_mm512_add_epi32(vblock, _mm512_set1_epi32(i+6)), &(text_int[0]), 4);
      vp7 = _mm512_mullo_epi32(vp7, vsize);
      vp7 = _mm512_add_epi32(vss, vp7);
      vss = _mm512_i32gather_epi32(vp7, trans_table, 4);

      __m512i vp8 = _mm512_i32gather_epi32(_mm512_add_epi32(vblock, _mm512_set1_epi32(i+7)), &(text_int[0]), 4);
      vp8 = _mm512_mullo_epi32(vp8, vsize);
      vp8 = _mm512_add_epi32(vss, vp8);
      vss = _mm512_i32gather_epi32(vp8, trans_table, 4);
    }

   ppp+=16;
    __mmask16 m = _mm512_cmpneq_epi32_mask(vss1, vss);
    if(m==0) {
      s = _mm512_mask_reduce_add_epi32(0x8000, vss);
      cur+=16*BLOCK;
      qqq+=15;
    } else {
      s = _mm512_mask_reduce_add_epi32(least_one[m], vss);
      cur += least_one_pos[m] * BLOCK;
      qqq+=least_one_pos[m];
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
