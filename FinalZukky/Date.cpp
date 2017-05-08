#include <iostream>
#include <algorithm>
using namespace std;


const int n_day[13] = {0, 31, 28, 31, 30, 31, 30, 31, 31, 30, 31, 30, 31};
 
struct Date {
  // year, month[1,12], day(1 origin), a day of the week[0,6]
  int y, m, d, e;
  Date() {}
  Date(int y, int m, int d, int e = -1) : y(y), m(m), d(d), e(e) {}
 
  bool operator < (const Date &o) const {
    if(y != o.y) return y < o.y;
    if(m != o.m) return m < o.m;
    return d < o.d;
  }
 
  bool isLeapYear() {
    if(y % 400 == 0) return true;
    if(y % 100 == 0) return false;
    if(y % 4 == 0) return true;
    return false;
  }

  int getDays() {
    int nd = n_day[m];
    if(m == 2) nd += isLeapYear();
    return nd;
  }

  void next() {
    e = (e+1)%7;
    if(d == getDays()) {
      d = 1;
      if(m == 12) {
        m = 1;
        ++y;
      } else {
        ++m;
      }
    } else {
      ++d;
    }
  }
 
  void back() {
    e = (e+6)%7;
    if(d == 1) {
      if(m == 1) {
        m = 12;
        --y;
      } else {
        --m;
      }
      d = getDays();
    } else {
      --d;
    }
  }
};

// AOJ 0125
int main() {
  int y1, m1, d1, y2, m2, d2;
  while(cin >> y1 >> m1 >> d1 >> y2 >> m2 >> d2) {
    if(y1 < 0 || m1 < 0 || d1 < 0 || y2 < 0 || m2 < 0 || d2 < 0)
      break;
    Date d(y1,m1,d1), end(y2,m2,d2);
    int res = 0;
    while(d < end) {
      ++res;
      d.next();
    }
    cout << res << endl;
  }
  return 0;
}
