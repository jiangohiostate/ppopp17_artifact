#include <iostream>
#include <fstream>
#include <sstream>
#include <map>
#include <sys/time.h>
using namespace std;

int num_states;
int num_alphs;
int *trans_table;
map<char, int> lookup_table;
int cur_in = 0;

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
  struct timeval tv1, tv2;
  struct timezone tz1, tz2;
  gettimeofday(&tv1, &tz1);
  for(int i=0;i<text.length();i++) {
    int index = text_int[i] * num_states + s;
    s = trans_table[index];
  }
  gettimeofday(&tv2, &tz2);
  cout << "time used: " << tv2.tv_usec - tv1.tv_usec + 1000000*(tv2.tv_sec - tv1.tv_sec) << "ms" << endl;
  cout << "final state: " << s << endl;
}


int main(int argc, char *argv[])
{
  input_rex(argv[1]);
  match(argv[2]);


  return 0;
}
