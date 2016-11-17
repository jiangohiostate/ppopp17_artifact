#include <iostream>
#include <string>
#include <fstream>
#include <vector>
#include <sstream>
#include <sys/time.h>
#include "immintrin.h"
using namespace std;

#define BLOCK 64

int num_states;
int *trans_table;

__m512i *spec_states;
__m512i *vind;

__mmask16 least_one[65536];
int least_one_pos[65536];

void print_vector(__m512i v) {
  for(int i=0;i<16;i++) {
    cout << *((int *)&v + i) << " ";
  }
  cout << endl;
}

struct HtmlToken {
  string token_str;
};

vector<HtmlToken> tokenize(string text)
{
  vector<HtmlToken> tokens;
  HtmlToken t;
  int s = 0;

  int cur = 0;
  float ppp = 0;
  float qqq = 0;

  __m512i vp;
  __m512i vs;
  __m512i vss;
  __m512i vtwo = _mm512_set1_epi32(2);
  int *vss_buffer = (int *)_mm_malloc(sizeof(int)*16, 64);

  int *text_int = (int *)_mm_malloc(sizeof(int)*text.length(), 64);
  for(int i=0;i<text.length();i++) {
    text_int[i] = (int)text[i] * num_states;
  }


  struct timeval tv1, tv2;
  struct timezone tz1, tz2;
  gettimeofday(&tv1, &tz1);

  while(cur + 8*BLOCK < text.length()) {
    vs = spec_states[s];
    vss = vs;

    for(int i=cur;i<cur+BLOCK;i++) {
      vp = _mm512_i32gather_epi32(vind[i-cur], text_int+cur, 4);
      vp = _mm512_add_epi32(vss, vp);
      vss = _mm512_i32gather_epi32(vp, trans_table, 4);
    }

    ppp+=7;
    _mm512_store_epi32(vss_buffer, vss);
    s = vss_buffer[0];
    if(s!=0) {
      cur+=BLOCK;
    } else {
      qqq++;
      s = vss_buffer[1];
      if(s!=0) {
        cur+=2*BLOCK;
      } else {
      qqq++;
        s = vss_buffer[4];
        if(s!=0) {
          cur+=3*BLOCK;
        } else {
      qqq++;
          s = vss_buffer[6];
          if(s!=0) {
            cur+=4*BLOCK;
          } else {
      qqq++;
            s = vss_buffer[8];
            if(s!=0) {
              cur+=5*BLOCK;
            } else {
      qqq++;
              s = vss_buffer[10];
              if(s!=0) {
                cur+=6*BLOCK;
              } else {
      qqq++;
                s = vss_buffer[12];
                if(s!=0) {
                  cur+=7*BLOCK;
                } else {
      qqq++;
                  s = vss_buffer[14];
                  cur+=8*BLOCK;
                }
              }
            }
          }
        }
      }
    }



  }


  for(int i=cur;i<text.length();i++) {
    s = trans_table[(text_int[i]) + s];
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
  cout << "used time: " <<  tv2.tv_usec - tv1.tv_usec + 1000000*(tv2.tv_sec - tv1.tv_sec) << endl;
  cout << "spec succ rate: " << (qqq+ppp)/(2*ppp) << endl;

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

  spec_states = (__m512i *)_mm_malloc(sizeof(__m512i)*num_states, 64);
  int *tmp_states = (int *)_mm_malloc(sizeof(int)*16, 64);
  for(int i=0;i<num_states;i++) {
    tmp_states[0] = i;
    for(int j=1;j<4;j++)tmp_states[j] = j-1;
    for(int j=4;j<16;j++)tmp_states[j] = j%2;
    spec_states[i] = _mm512_load_epi32(tmp_states);
  }

  vind = (__m512i *)_mm_malloc(sizeof(__m512i)*BLOCK, 64);
  for(int i=0;i<BLOCK;i++) {
    tmp_states[0] = i;
    for(int j=1;j<4;j++)tmp_states[j] = i+BLOCK;
    for(int j=4;j<16;j++)tmp_states[j] = i+j/2*BLOCK;
    vind[i] = _mm512_load_epi32(tmp_states);
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
