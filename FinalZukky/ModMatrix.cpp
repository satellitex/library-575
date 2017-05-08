//validate by AOJ 2397
#include <bits/stdc++.h>
using namespace std;
typedef long long lli;
typedef lli elm;
typedef vector<vector<elm> > mat;

mat getE(int n) {
  mat e(n, vector<elm>(n, 0));
  for(int i = 0; i < n; ++i) {
    e[i][i] = 1;
  }
  return e;
}

mat mod_multiply(const mat &a, const mat &b, const elm &mod) {
  if(a.size() == 0 || b.size() == 0 ||
     a[0].size() != b.size()) return mat();
  int s = a.size();
  int t = b.size();
  int u = b[0].size();
  mat c(s, vector<elm>(u, 0));
  for(int i = 0; i < s; ++i) {
    for(int j = 0; j < u; ++j) {
      for(int k = 0; k < t; ++k) {
        ( c[i][j] += a[i][k] * b[k][j] % mod ) % mod;
      }
    }
  }
  return c;
}

mat mod_pow(const mat &a, const lli &n, const elm &mod) {
  if(n == 0) return getE(a.size());
  mat b = mod_pow(mod_multiply(a, a, mod), n/2, mod);
  if(n & 1) return mod_multiply(b, a, mod);
  else return b;
}

// by Spaghetti Source
// O( n^3 ).
elm mod_det(mat a, const elm &mod) {
  const int n = a.size();

  for(int i = 0; i < n; ++i) {
    for(int j = 0; j < n; ++j) {
      a[i][j] = a[i][j] % mod;
      if(a[i][j] < 0) a[i][j] += mod;
    }
  }

  elm res = 1;
  for (int i = 0; i < n; ++i) {
    int pivot = i;
    for (int j = i+1; j < n; ++j)
      if (abs(a[j][i]) > abs(a[pivot][i])) pivot = j;
    swap(a[pivot], a[i]);
    ( res *= a[i][i] * (i != pivot ? -1 : 1) ) %= mod;
    if(res < 0) res += mod;
    if (a[i][i] == 0) break;
    elm inv_a = mod_inv(a[i][i], mod);
    for(int j = i+1; j < n; ++j)
      for(int k = n-1; k >= i; --k) {
        ( a[j][k] -= a[i][k] * a[j][i] % mod * inv_a % mod ) %= mod;
        if(a[j][k] < 0) a[j][k] += mod;
      }
  }
  return res;
}

// by Spaghetti Source
// O( n^3 ).
/*
int rank(mat a) {
  const int n = a.size(), m = a[0].size();
  int r = 0;
  for (int i = 0; r < n && i < m; ++i) {
    int pivot = r;
    for (int j = r+1; j < n; ++j)
      if (abs(a[j][i]) > abs(a[pivot][i])) pivot = j;
    swap(a[pivot], a[r]);
    if (abs(a[r][i]) < eps) continue;
    for (int k = m-1; k >= i; --k)
      a[r][k] /= a[r][i];
    for(int j = r+1; j < n; ++j)
      for(int k = i; k < m; ++k)
        a[j][k] -= a[r][k] * a[j][i];
    ++r;
  }
  return r;
}
*/

int main() {

  return 0;
}
