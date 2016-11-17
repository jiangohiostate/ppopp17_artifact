#include <iostream>
#include <string>
#include <fstream>
#include <vector>
#include <sstream>
#include <sys/time.h>
using namespace std;

int num_states;
int *trans_table;


struct HtmlToken {
  string token_str;
};

vector<HtmlToken> tokenize(string text)
{
  vector<HtmlToken> tokens;
  HtmlToken t;
  int s = 0;

  struct timeval tv1, tv2;
  struct timezone tz1, tz2;
  gettimeofday(&tv1, &tz1);
  for(int i=0;i<text.length();i++) {
    int index = text[i] * num_states + s;
    s = trans_table[index];
   /* if(s==0)t.token_str = "";
    else {
      t.token_str += text[i];
      if(trans_table[127*num_states+s] < 0) {
        t.token_str += '>';
        tokens.push_back(t);
        s = 0;
      }
    }*/
  }
  gettimeofday(&tv2, &tz2);
  cout << "used time: " << tv2.tv_usec - tv1.tv_usec + 1000000*(tv2.tv_sec - tv1.tv_sec) << endl;

  cout << "final state: " << s << endl;

  return tokens;
}

void input_FSM(char *fsmfile) 
{
  ifstream fin(fsmfile);
  string line;
  getline(fin, line);
  stringstream sin1(line);
  sin1 >> num_states;

  trans_table = (int *)_mm_malloc(sizeof(int)*num_states*128, 64);
  for(int i=0;i<128;i++) {
    for(int j=0;j<num_states;j++) {
      trans_table[i*num_states+j] = 0;
    }
  }

  // the last row of the states stores the types of the token
  getline(fin, line);
  stringstream sin2(line);
  int fstate;
  while(sin2 >> fstate) {
    trans_table[127*num_states+fstate] = -1 * fstate;
  }


  int s, e;
  char c;

  while(getline(fin, line)) {
    stringstream sin(line);
    sin >> s >> c >> e;
    trans_table[c*num_states+s] = e;
  }

}


int main(int argc, char *argv[])
{
  ifstream fin(argv[1]); 
  string line, html;
  while(getline(fin, line)) {
    for(int i=0;i<line.length();i++) {
      if(line[i] >='a' && line[i]<='z') line[i] = line[i] - 'a' + 'A';
      if(line[i] >= 127 || line[i] < 0)line[i]=0;
    }
    html += line;
  }
  input_FSM("datasets/html.fsm");
  vector<HtmlToken> tokens = tokenize(html);
  for(auto t : tokens) {
    cout << t.token_str << endl;
  }

  return 0;
}
