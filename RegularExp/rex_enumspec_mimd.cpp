#include <iostream>
#include <fstream>
#include <sstream>
#include <map>
#include <sys/time.h>
#include <stdlib.h>
#include "immintrin.h"
#include <limits.h>

using namespace std;

int num_states;
int num_alphs;
int *trans_table;
map<char, int> lookup_table;
int cur_in = 0;
int *spec_table_temp;

__mmask16 least_one[65536];
int least_one_pos[65536];

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
 //init the pos of the first 1 from least significant bit
  for(int i=1;i<65536;i++) {
    int c = 0;
    int t = i;
    while(t % 2 == 0) {
      t /= 2;
      c++;
    }
    least_one[i] = (1 << c);
    least_one_pos[i] = c;
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
 spec_table_temp = (int *)_mm_malloc(sizeof(int)*16*(num_alphs+1)*(num_alphs+1), 64);

 int *all_states = (int *)_mm_malloc(sizeof(int)*num_states, 64);
 int *to_states = (int *)_mm_malloc(sizeof(int)*num_states, 64);
 int *scount = new int[num_states];
 int max, p;


 for(int i1=0;i1<=num_alphs;i1++) {
   for(int i2=0;i2<=num_alphs;i2++) {
     for(int k=0;k<num_states;k++) {
       to_states[k] = trans_table[trans_table[i1*num_states+k]+i2*num_states]; 
     }
     for(int k=0;k<num_states;k++) {
       scount[k] = 0; 
     }

     for(int k=0;k<num_states;k++) {
       scount[to_states[k]]++; 
     }

     for(int j=0;j<16;j++) {
       max = INT_MIN;
       for(int k=0;k<num_states;k++) {
         if(scount[k] > max) {
           max = scount[k];
           p = k;
         }
       }
       scount[p] = -1;
       spec_table_temp[16*(i1*(num_alphs+1)+i2)+j] = p;
     }
   }
 }


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

  int num_blocks = text.length() / BLOCK;
  __m512i *res_states = (__m512i *) _mm_malloc(sizeof(__m512i)*num_blocks, 64);
  __m512i *spec_states = (__m512i *) _mm_malloc(sizeof(__m512i)*num_blocks, 64);
  __m512i vperm = _mm512_set_epi32(7,6,5,4,3,2,1,0,15,14,13,12,11,10,9,8);

  struct timezone tz1, tz2;
  struct timeval tv1, tv2;
  gettimeofday(&tv1, &tz1);

  int submitted = 0;

  int nb = num_blocks % 2 == 0 ? num_blocks-1:num_blocks;
#pragma omp parallel for num_threads(nth) 
  for(int j = 1; j < nb; j+=2) {
    int start = j * BLOCK;
    int end = start + BLOCK;
    int index1 = text_int[start-2]*(num_alphs+1)+text_int[start-1];
    int index2 = text_int[end-2]*(num_alphs+1)+text_int[end-1];
    __m512i ss1 = _mm512_load_epi32(spec_table_temp+index1*16);
    __m512i ss2 = _mm512_load_epi32(spec_table_temp+index2*16);

    spec_states[j] = ss1;
    spec_states[j+1] = ss2;
    __m512i vs = _mm512_mask_permutevar_epi32(ss1, 0xFF00, vperm, ss2);

    for(int i=start; i<end; i+=8) {
      __m512i a1 = _mm512_set1_epi32(text_int[i]);
      __m512i b1 = _mm512_set1_epi32(text_int[i+BLOCK]);
      __m512i vi1 = _mm512_mask_permutevar_epi32(a1, 0xFF00, vperm, b1);
      vi1 = _mm512_mullo_epi32(vi1, _mm512_set1_epi32(num_states));
      vi1 = _mm512_add_epi32(vs, vi1);
      vs = _mm512_i32gather_epi32(vi1, trans_table, 4);

      __m512i a2 = _mm512_set1_epi32(text_int[i+1]);
      __m512i b2 = _mm512_set1_epi32(text_int[i+1+BLOCK]);
      __m512i vi2 = _mm512_mask_permutevar_epi32(a2, 0xFF00, vperm, b2);
      vi2 = _mm512_mullo_epi32(vi2, _mm512_set1_epi32(num_states));
      vi2 = _mm512_add_epi32(vs, vi2);
      vs = _mm512_i32gather_epi32(vi2, trans_table, 4);

      __m512i a3 = _mm512_set1_epi32(text_int[i+2]);
      __m512i b3 = _mm512_set1_epi32(text_int[i+2+BLOCK]);
      __m512i vi3 = _mm512_mask_permutevar_epi32(a3, 0xFF00, vperm, b3);
      vi3 = _mm512_mullo_epi32(vi3, _mm512_set1_epi32(num_states));
      vi3 = _mm512_add_epi32(vs, vi3);
      vs = _mm512_i32gather_epi32(vi3, trans_table, 4);

      __m512i a4 = _mm512_set1_epi32(text_int[i+3]);
      __m512i b4 = _mm512_set1_epi32(text_int[i+3+BLOCK]);
      __m512i vi4 = _mm512_mask_permutevar_epi32(a4, 0xFF00, vperm, b4);
      vi4 = _mm512_mullo_epi32(vi4, _mm512_set1_epi32(num_states));
      vi4 = _mm512_add_epi32(vs, vi4);
      vs = _mm512_i32gather_epi32(vi4, trans_table, 4);

      __m512i a5 = _mm512_set1_epi32(text_int[i+4]);
      __m512i b5 = _mm512_set1_epi32(text_int[i+4+BLOCK]);
      __m512i vi5 = _mm512_mask_permutevar_epi32(a5, 0xFF00, vperm, b5);
      vi5 = _mm512_mullo_epi32(vi5, _mm512_set1_epi32(num_states));
      vi5 = _mm512_add_epi32(vs, vi5);
      vs = _mm512_i32gather_epi32(vi5, trans_table, 4);

      __m512i a6 = _mm512_set1_epi32(text_int[i+5]);
      __m512i b6 = _mm512_set1_epi32(text_int[i+5+BLOCK]);
      __m512i vi6 = _mm512_mask_permutevar_epi32(a6, 0xFF00, vperm, b6);
      vi6 = _mm512_mullo_epi32(vi6, _mm512_set1_epi32(num_states));
      vi6 = _mm512_add_epi32(vs, vi6);
      vs = _mm512_i32gather_epi32(vi6, trans_table, 4);

      __m512i a7 = _mm512_set1_epi32(text_int[i+6]);
      __m512i b7 = _mm512_set1_epi32(text_int[i+6+BLOCK]);
      __m512i vi7 = _mm512_mask_permutevar_epi32(a7, 0xFF00, vperm, b7);
      vi7 = _mm512_mullo_epi32(vi7, _mm512_set1_epi32(num_states));
      vi7 = _mm512_add_epi32(vs, vi7);
      vs = _mm512_i32gather_epi32(vi7, trans_table, 4);

      __m512i a8 = _mm512_set1_epi32(text_int[i+7]);
      __m512i b8 = _mm512_set1_epi32(text_int[i+7+BLOCK]);
      __m512i vi8 = _mm512_mask_permutevar_epi32(a8, 0xFF00, vperm, b8);
      vi8 = _mm512_mullo_epi32(vi8, _mm512_set1_epi32(num_states));
      vi8 = _mm512_add_epi32(vs, vi8);
      vs = _mm512_i32gather_epi32(vi8, trans_table, 4);



    }
    res_states[j] = vs;
    res_states[j+1] = _mm512_permutevar_epi32(vperm, vs);
  }

  if(num_blocks % 2 == 0) {
    int start = (num_blocks-1) * BLOCK;
    int end = start + BLOCK;
    int index = text_int[start-2]*(num_alphs+1)+text_int[start-1];
    __m512i vs = _mm512_load_epi32(spec_table_temp+index*16);

    spec_states[num_blocks-1] = vs;

    for(int i=start; i<end; i++) {
      __m512i vi = _mm512_set1_epi32(text_int[i]);
      vi = _mm512_mullo_epi32(vi, _mm512_set1_epi32(num_states));
      vi = _mm512_add_epi32(vs, vi);
      vs = _mm512_i32gather_epi32(vi, trans_table, 4);
    }
    res_states[num_blocks-1] = vs;
  }

  int cur = 0;
  int ppp = 0;

  while(cur < text.length()) {
    int end = cur + BLOCK < text.length() ? cur + BLOCK : text.length();
    for(int i=cur;i<end;i++) {
      s = trans_table[s+(text_int[i])*num_states];
    }
    ppp++;
    submitted++;
    cur += BLOCK;
    while(submitted < num_blocks) {
      __m512i vt = _mm512_set1_epi32(s);
      __mmask16 m = _mm512_mask_cmpeq_epi32_mask(0x00FF, vt, spec_states[submitted]);
      if(m) {
        s = _mm512_mask_reduce_add_epi32(least_one[m], res_states[submitted]); 
        submitted++;
        cur += BLOCK;
      }
      else {
        break; 
      }
    }
  }

  gettimeofday(&tv2, &tz2);
  cout << "final state: " << s << endl;
  cout << "time used: " << tv2.tv_usec - tv1.tv_usec + 1000000 * (tv2.tv_sec - tv1.tv_sec) << "ms"<<endl;
  cout << "spec rate: " << 1.0*(num_blocks - ppp) / (num_blocks - 1) << endl;

}


int main(int argc, char *argv[])
{
  int nth = atoi(argv[3]);
  input_rex(argv[1]);
  match(argv[2], nth);


  return 0;
}
