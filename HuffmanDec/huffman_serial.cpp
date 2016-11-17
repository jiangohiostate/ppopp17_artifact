#include <iostream>
#include <map>
#include <fstream>
#include <climits>
#include <vector>
#include <iomanip>
#include <bitset>
#include <sys/time.h>
using namespace std;

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

struct StateTuple {
  int s[2];
  int c = EOF;
};

string compressed;
HuffmanNode *tree;
int* stable;
int nstates;


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
  nstates = 2*size-1;
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
  store_code(tree);

  fin.close();

  int osize = 0;
  ifstream fin2(filename);
  while((c = fin2.get())!=EOF) {
    compressed += hist[c]->code;
    osize += 8;
  }
  
  cout << "original_size: " << osize << " bits" << endl;
  cout << "huffman_size: " << compressed.length() << " bits" << endl;
}

void decode_serial()
{
  int s = 0;
  cout << "number of states: " << nstates << endl;
  struct timezone tz1, tz2;
  struct timeval tv1, tv2;
  gettimeofday(&tv1, &tz1);

  int i=0;
  for(;i<compressed.length();i++) {
    s = stable[s+(compressed[i]-'0')*nstates];
  } 
  gettimeofday(&tv2, &tz2);
  cout << "final state: " << s << endl;
 cout << "time used: " << tv2.tv_usec - tv1.tv_usec + 1000000 * (tv2.tv_sec - tv1.tv_sec) << "ms"<< endl;
}

int main(int argc, char* argv[])
{
  encode(argv[1]);
  decode_serial();
  return 0;
}
