#include <bits/stdc++.h>
using namespace std;

typedef long long lli;

// O(log max(a,b))
lli extgcd(lli a, lli b, lli &x, lli &y) {
  lli d = a;
  if(b != 0) {
    d = extgcd(b, a % b, y, x);
    y -= (a / b) * x;
  } else {
    x = 1; y = 0;
  }
  return d;
}

lli mod_inverse(lli a, lli m) {
  lli x, y;
  extgcd(a, m, x, y);
  return (m + x % m) % m;
}
