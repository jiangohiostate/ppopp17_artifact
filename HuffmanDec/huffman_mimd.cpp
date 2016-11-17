#include <iostream>
#include <map>
#include <fstream>
#include <climits>
#include <vector>
#include <iomanip>
#include <bitset>
#include <sys/time.h>
#include "immintrin.h"
#include <omp.h>
using namespace std;

#define BLOCK 16384

struct HuffmanNode {
  HuffmanNode(char _c, int _count, struct HuffmanNode* p1, struct HuffmanNode *p2):c(_c), count(_count), left(p1), right(p2){}
  char c;
  int id;
  int count;
  struct HuffmanNode* left;
  struct HuffmanNode* right;
  string code;
};

void print_vec(__m512i v)
{
  for(int i=0;i<16;i++) {
    cout << *((int *)&v + i) << " ";
  }
  cout << endl;
}


void store_code(struct HuffmanNode* p, string prefix="")
{
  if(p != NULL) {
    if(p->left == NULL && p->right == NULL)p->code = prefix;
    if(p->left) store_code(p->left, prefix+"0");
    if(p->right) store_code(p->right, prefix+"1");
  }
}

string compressed;
vector<int> compressed_int;
HuffmanNode *tree;
int* stable;
__m512i *spec_mm_table;
int *input_table;
__m512i input_mm_table[4];
int size_of_states;
__mmask16 least_one[65536];
int least_one_pos[65536];
__m512i vtsize;
int *spec_table_temp;


void encode(string filename) {
  // count number of each character
  ifstream fin(filename.c_str());
  map<char, HuffmanNode*> hist;
  char c;
  while((c=fin.get())!=EOF) {
    map<char, HuffmanNode*>::iterator it = hist.find(c);
    if(it!=hist.end()) {
      (it->second->count)++;
    } else {
      hist[c] = new HuffmanNode(c, 1, NULL, NULL);
    }
  }
   // store the nodes in a vector
  vector<HuffmanNode*> vhist;
  for(auto &kv : hist) {
    vhist.push_back(kv.second);
  }
  // build the huffman tree
  int size = hist.size();
  int min = INT_MAX;
  vector<HuffmanNode*>::iterator itt;
  vector<HuffmanNode*> state_list;
  HuffmanNode *n3;
  for(int i=0;i<size-1;i++) {
    min = INT_MAX;
    for(vector<HuffmanNode*>::iterator it=vhist.begin();it!=vhist.end();it++) {
      if((*it)->count < min) {
        min = (*it)->count;
        itt = it;
      }
    }
    HuffmanNode *n1 = *itt;
    state_list.push_back(n1);
    vhist.erase(itt);
    min = INT_MAX;
    for(vector<HuffmanNode*>::iterator it=vhist.begin();it!=vhist.end();it++) {
      if((*it)->count < min) {
        min = (*it)->count;
        itt = it;
      }
    }
    HuffmanNode *n2 = *itt;
    state_list.push_back(n2);
    vhist.erase(itt);
    n3 = new HuffmanNode(NULL, n1->count+n2->count, n1, n2);
    vhist.push_back(n3);
  }
  state_list.push_back(n3);
  tree = n3;


  //build state transfer table
  size_of_states = 2*size-1;
  vtsize = _mm512_set1_epi32(2*size-1);
  stable = (int *)_mm_malloc(sizeof(int)*3*(2*size-1), 64);
  int cur_state = 0;
  for(int i=state_list.size()-1;i>=0;i--) {
    state_list[i]->id = cur_state++;
  }
  for(int i=state_list.size()-1;i>=0;i--) {
    if(state_list[i]->left != NULL) {
      //non_leaf
      stable[(state_list[i]->id)] = state_list[i]->left->id;
      stable[(state_list[i]->id)+2*size-1] = state_list[i]->right->id;
    } else {
      //leaf
      stable[(state_list[i]->id)] = stable[0];
      stable[(state_list[i]->id)+2*size-1] = stable[2*size-1];
      stable[(state_list[i]->id)+4*size-2] = state_list[i]->c;
    }
  }
  int max, p;
  int* temp0 = new int[2*size-1];
  int* temp1 = new int[2*size-1];
  for(int i=0;i<2*size-1;i++) {
   temp0[i] = stable[i];  
   temp1[i] = stable[i+2*size-1];  
  }

  spec_table_temp = (int *)_mm_malloc(sizeof(int)*4*16, 64);
  int *count = new int[2*size-1];
  for(int i=0;i<2*size-1;i++) {
   count[i] = 0; 
  }
  for(int i=0;i<2*size-1;i++) {
   count[temp0[temp0[i]]]++; 
  }
  for(int t=0;t<16;t++) {
    max = INT_MIN;
    for(int i=0;i<2*size-1;i++) {
      if(count[i] > max) {
        max = count[i];
        p = i;
      }
    }
    spec_table_temp[t] = p;
    count[p] = -1;
  }

  for(int i=0;i<2*size-1;i++) {
    count[i] = 0; 
  }
  for(int i=0;i<2*size-1;i++) {
    count[temp1[temp0[i]]]++; 
  }
  for(int t=0;t<16;t++) {
    max = INT_MIN;
    for(int i=0;i<2*size-1;i++) {
      if(count[i] > max) {
        max = count[i];
        p = i;
      }
    }
    spec_table_temp[16+t] = p;
    count[p] = -1;
  }
  for(int i=0;i<2*size-1;i++) {
    count[i] = 0; 
  }
  for(int i=0;i<2*size-1;i++) {
    count[temp0[temp1[i]]]++; 
  }
  for(int t=0;t<16;t++) {
    max = INT_MIN;
    for(int i=0;i<2*size-1;i++) {
      if(count[i] > max) {
        max = count[i];
        p = i;
      }
    }
    spec_table_temp[32+t] = p;
    count[p] = -1;
  }
  for(int i=0;i<2*size-1;i++) {
    count[i] = 0; 
  }
  for(int i=0;i<2*size-1;i++) {
    count[temp1[temp1[i]]]++; 
  }
  for(int t=0;t<16;t++) {
    max = INT_MIN;
    for(int i=0;i<2*size-1;i++) {
      if(count[i] > max) {
        max = count[i];
        p = i;
      }
    }
    spec_table_temp[48+t] = p;
    count[p] = -1;
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

  store_code(tree);

  fin.close();

  int osize = 0;
  ifstream fin2(filename);
  while((c = fin2.get())!=EOF) {
    compressed += hist[c]->code;
    osize += 8;
  }
  for(int i=0;i<compressed.length();i++) {
    compressed_int.push_back(compressed[i]-'0');
  }

  cout << "original_size: " << osize << " bits" << endl;
  cout << "huffman_size: " << compressed.length() << " bits" << endl;
}

__m512i* spec_states;
__m512i* res_states;

void decode_mimd(int nth)
{
  int num_blocks = compressed_int.size() / BLOCK;
  spec_states = (__m512i *) _mm_malloc(sizeof(__m512i)*num_blocks, 64);
  res_states = (__m512i *) _mm_malloc(sizeof(__m512i)*num_blocks, 64);

  struct timezone tz1, tz2;
  struct timeval tv1, tv2;
  gettimeofday(&tv1, &tz1);

  int submitted = 0;

#pragma omp parallel for num_threads(nth) 
  for(int j = 1; j < num_blocks; j++) {
    int start = j * BLOCK;
    int end = start + BLOCK;
    int t = (compressed_int[start - 2]) * 2 + compressed_int[start - 1];
    __m512i vs = _mm512_load_epi32(t*16 + spec_table_temp);
    spec_states[j] = vs;
    for(int i=start; i<end; i+=8) {
      __m512i vi1 = _mm512_set1_epi32(compressed_int[i]);
      vi1 = _mm512_mullo_epi32(vi1, _mm512_set1_epi32(size_of_states));
      vi1 = _mm512_add_epi32(vs, vi1);
      vs = _mm512_i32gather_epi32(vi1, stable, 4);

      __m512i vi2 = _mm512_set1_epi32(compressed_int[i+1]);
      vi2 = _mm512_mullo_epi32(vi2, _mm512_set1_epi32(size_of_states));
      vi2 = _mm512_add_epi32(vs, vi2);
      vs = _mm512_i32gather_epi32(vi2, stable, 4);

      __m512i vi3 = _mm512_set1_epi32(compressed_int[i+2]);
      vi3 = _mm512_mullo_epi32(vi3, _mm512_set1_epi32(size_of_states));
      vi3 = _mm512_add_epi32(vs, vi3);
      vs = _mm512_i32gather_epi32(vi3, stable, 4);

      __m512i vi4 = _mm512_set1_epi32(compressed_int[i+3]);
      vi4 = _mm512_mullo_epi32(vi4, _mm512_set1_epi32(size_of_states));
      vi4 = _mm512_add_epi32(vs, vi4);
      vs = _mm512_i32gather_epi32(vi4, stable, 4);

      __m512i vi5 = _mm512_set1_epi32(compressed_int[i+4]);
      vi5 = _mm512_mullo_epi32(vi5, _mm512_set1_epi32(size_of_states));
      vi5 = _mm512_add_epi32(vs, vi5);
      vs = _mm512_i32gather_epi32(vi5, stable, 4);

      __m512i vi6 = _mm512_set1_epi32(compressed_int[i+5]);
      vi6 = _mm512_mullo_epi32(vi6, _mm512_set1_epi32(size_of_states));
      vi6 = _mm512_add_epi32(vs, vi6);
      vs = _mm512_i32gather_epi32(vi6, stable, 4);

      __m512i vi7 = _mm512_set1_epi32(compressed_int[i+6]);
      vi7 = _mm512_mullo_epi32(vi7, _mm512_set1_epi32(size_of_states));
      vi7 = _mm512_add_epi32(vs, vi7);
      vs = _mm512_i32gather_epi32(vi7, stable, 4);

      __m512i vi8 = _mm512_set1_epi32(compressed_int[i+7]);
      vi8 = _mm512_mullo_epi32(vi8, _mm512_set1_epi32(size_of_states));
      vi8 = _mm512_add_epi32(vs, vi8);
      vs = _mm512_i32gather_epi32(vi8, stable, 4);





    }
    res_states[j] = vs;
  }

 
  int s = 0;
  int cur = 0;
  int ppp = 0;

  while(cur < compressed_int.size()) {
    int end = cur + BLOCK < compressed_int.size() ? cur + BLOCK : compressed_int.size();
    for(int i=cur;i<end;i++) {
      s = stable[s+(compressed_int[i])*size_of_states];
    }
    ppp++;
    submitted++;
    cur += BLOCK;
    while(submitted < num_blocks) {
      __m512i vt = _mm512_set1_epi32(s);
      __mmask16 m = _mm512_cmpeq_epi32_mask(vt, spec_states[submitted]);
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
  cout << "time used: " << tv2.tv_usec - tv1.tv_usec + 1000000 * (tv2.tv_sec - tv1.tv_sec) << endl;
  cout << "spec rate: " << 1.0*(num_blocks - ppp) / (num_blocks - 1) << endl;
}


int main(int argc, char* argv[])
{
  encode(argv[1]);
  //cout << compressed << endl;
  int nth = atoi(argv[2]);
  decode_mimd(nth);
  return 0;
}
