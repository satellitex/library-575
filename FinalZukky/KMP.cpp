#include <iostream>
#include <algorithm>
#include <cctype>
#include <string>
#include <vector>
using namespace std;

class KMP {
  string pattern;
  vector<int> next;
public:
  KMP(string p) : pattern(p), next(vector<int>(p.size()+1)) {
    next[0] = -1;
    for(int i = 0, j = -1; i < pattern.size(); ++i, ++j, next[i] = j)
      while(j >= 0 && pattern[i] != pattern[j]) j = next[j];
  }

  int search(const string &text) {
    int i, j;
    i = j = 0;
    while(i < text.size() && j < pattern.size()) {
      if(text[i] == pattern[j]) {
        ++i;
        ++j;
      } else if(j == 0) {
        ++i;
        j = 0;
      } else {
        j = next[j];
      }
    }
    if(j == pattern.size()) return i-pattern.size();
    else return -1;
  }

  vector<int> searchV(const string &text) {
    vector<int> res;
    int i, j;
    i = j = 0;
    while(i < text.size() && j < pattern.size()) {
      if(text[i] == pattern[j]) {
        ++i;
        ++j;
      } else if(j == 0){
        ++i;
        j = 0;
      } else {
        j = next[j];
      }
    }
    if(j == pattern.size()){
      res.push_back(i - j);
      i = i - j + 1;
      j = 0;
    }
    return res;
  }

  int search2(const string &text) { // uva 11475
    int i, j;
    i = j = 0;
    while(i < text.size() && j < pattern.size()) {
      if(text[i] == pattern[j]) {
        ++i;
        ++j;
      } else if(j == 0) {
        ++i;
        j = 0;
      } else {
        j = next[j];
      }
    }
    return j;
  }
};

// uva 11475
int main() {
  string str;
  while(cin >> str) {
    string rstr = str;
    reverse(rstr.begin(), rstr.end());
    KMP kmp(rstr);
    cout << str << rstr.substr(kmp.search2(str)) << endl;
  }
  return 0;
}
