#include <iostream>
#include <fstream>
#include <sstream>
#include <string>
#include <sys/time.h>
using namespace std;

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
}


bool check_div(string num) {
  int *num_int = (int *)_mm_malloc(sizeof(int)*num.length(), 64);
  for(int i=0;i<num.length();i++) {
    num_int[i] = num[i] - '0';
  }

  int s = 0;
  struct timeval tv1, tv2;
  struct timezone tz1, tz2;
  gettimeofday(&tv1, &tz1);
  for(int i=0;i<num.length();i++) {
    int index = (num_int[i]) * num_states + s;
    s = trans_table[index];
  }
  gettimeofday(&tv2, &tz2);
  cout << "used time: " << tv2.tv_usec - tv1.tv_usec + 1000000 * (tv2.tv_sec - tv1.tv_sec) << endl;
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
