#include <iostream>
#include <map>
#include <fstream>
#include <climits>
#include <vector>
#include <iomanip>
#include <bitset>
#include <sys/time.h>
#include "immintrin.h"
using namespace std;

#define BLOCK 4096

struct HuffmanNode {
  HuffmanNode(char _c, int _count, struct HuffmanNode* p1, struct HuffmanNode *p2):c(_c), count(_count), left(p1), right(p2){}
  char c;
  int id;
  int count;
  struct HuffmanNode* left;
  struct HuffmanNode* right;
  string code;
};

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
int *spec_table;
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
  for(int i=2;i<65536;i++) {
    int c = 0;
    int t = i & 0xFFFE;
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


void decode_simd()
{
  int cur = 0;
  int s = 0;
  __m512i vs;
  __m512i vss;
  __m512i vone = _mm512_set1_epi32(0);
  __m512i vsize = _mm512_set1_epi32(size_of_states);

  __m512i vperm = _mm512_set_epi32(7,6,5,4,3,2,1,0,15,14,13,12,11,10,9,8);
  __m512i vperm1 = _mm512_set_epi32(14,13,12,11,10,9,8,7,6,5,4,3,2,1,0,15);

  float ppp = 0, qqq = 0;
  struct timezone tz1, tz2;
  struct timeval tv1, tv2;
  gettimeofday(&tv1, &tz1);

  while(cur + 2 * BLOCK < compressed.length()) {
    int index = (compressed_int[cur+BLOCK-2])*2+compressed_int[cur+BLOCK-1];
    __m512i vs1 = _mm512_load_epi32(spec_table_temp+16*index);
    vs = _mm512_set1_epi32(s);
    vs = _mm512_mask_permutevar_epi32(vs, 0xFFFE, vperm1, vs1);
    vss = vs;
    for(int i=cur;i<cur+BLOCK;i+=8) {
      __m512i a1 = _mm512_set1_epi32(compressed_int[i]);
      __m512i b1 = _mm512_set1_epi32(compressed_int[i+BLOCK]);
      __m512i vp1 = _mm512_mask_permutevar_epi32(a1, 0xFFFE, vperm, b1);
      vp1 = _mm512_mullo_epi32(vp1, vsize);
      vp1 = _mm512_add_epi32(vss, vp1);
      vss = _mm512_i32gather_epi32(vp1, stable, 4);

      __m512i a2 = _mm512_set1_epi32(compressed_int[i+1]);
      __m512i b2 = _mm512_set1_epi32(compressed_int[i+1+BLOCK]);
      __m512i vp2 = _mm512_mask_permutevar_epi32(a2, 0xFFFE, vperm, b2);
      vp2 = _mm512_mullo_epi32(vp2, vsize);
      vp2 = _mm512_add_epi32(vss, vp2);
      vss = _mm512_i32gather_epi32(vp2, stable, 4);

      __m512i a3 = _mm512_set1_epi32(compressed_int[i+2]);
      __m512i b3 = _mm512_set1_epi32(compressed_int[i+2+BLOCK]);
      __m512i vp3 = _mm512_mask_permutevar_epi32(a3, 0xFFFE, vperm, b3);
      vp3 = _mm512_mullo_epi32(vp3, vsize);
      vp3 = _mm512_add_epi32(vss, vp3);
      vss = _mm512_i32gather_epi32(vp3, stable, 4);

      __m512i a4 = _mm512_set1_epi32(compressed_int[i+3]);
      __m512i b4 = _mm512_set1_epi32(compressed_int[i+3+BLOCK]);
      __m512i vp4 = _mm512_mask_permutevar_epi32(a4, 0xFFFE, vperm, b4);
      vp4 = _mm512_mullo_epi32(vp4, vsize);
      vp4 = _mm512_add_epi32(vss, vp4);
      vss = _mm512_i32gather_epi32(vp4, stable, 4);

      __m512i a5 = _mm512_set1_epi32(compressed_int[i+4]);
      __m512i b5 = _mm512_set1_epi32(compressed_int[i+4+BLOCK]);
      __m512i vp5 = _mm512_mask_permutevar_epi32(a5, 0xFFFE, vperm, b5);
      vp5 = _mm512_mullo_epi32(vp5, vsize);
      vp5 = _mm512_add_epi32(vss, vp5);
      vss = _mm512_i32gather_epi32(vp5, stable, 4);

      __m512i a6 = _mm512_set1_epi32(compressed_int[i+5]);
      __m512i b6 = _mm512_set1_epi32(compressed_int[i+5+BLOCK]);
      __m512i vp6 = _mm512_mask_permutevar_epi32(a6, 0xFFFE, vperm, b6);
      vp6 = _mm512_mullo_epi32(vp6, vsize);
      vp6 = _mm512_add_epi32(vss, vp6);
      vss = _mm512_i32gather_epi32(vp6, stable, 4);

      __m512i a7 = _mm512_set1_epi32(compressed_int[i+6]);
      __m512i b7 = _mm512_set1_epi32(compressed_int[i+6+BLOCK]);
      __m512i vp7 = _mm512_mask_permutevar_epi32(a7, 0xFFFE, vperm, b7);
      vp7 = _mm512_mullo_epi32(vp7, vsize);
      vp7 = _mm512_add_epi32(vss, vp7);
      vss = _mm512_i32gather_epi32(vp7, stable, 4);

      __m512i a8 = _mm512_set1_epi32(compressed_int[i+7]);
      __m512i b8 = _mm512_set1_epi32(compressed_int[i+7+BLOCK]);
      __m512i vp8 = _mm512_mask_permutevar_epi32(a8, 0xFFFE, vperm, b8);
      vp8 = _mm512_mullo_epi32(vp8, vsize);
      vp8 = _mm512_add_epi32(vss, vp8);
      vss = _mm512_i32gather_epi32(vp8, stable, 4);
    }
    ppp+=2;
    qqq++;
    s =  _mm512_mask_reduce_add_epi32(0x0001, vss);
    int t = _mm512_mask_reduce_add_epi32(0x0002, vs);
    cur+=BLOCK;
    if(s==t) {
      qqq++;
      s = _mm512_mask_reduce_add_epi32(0x0002, vss);
      cur+=BLOCK;
    }
  }

  for(int i=cur;i<compressed_int.size();i++) {
    s = stable[s+(compressed_int[i])*size_of_states];
  }

  gettimeofday(&tv2, &tz2);
  cout << "final state: " << s << endl;
  cout << "time used: " << tv2.tv_usec - tv1.tv_usec + 1000000 * (tv2.tv_sec - tv1.tv_sec) << "ms" << endl;
  cout << "spec rate: " << qqq / ppp << endl;
}


int main(int argc, char* argv[])
{
  encode(argv[1]);
  decode_simd();
  return 0;
}
