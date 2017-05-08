//validate by AOJ 2397
#include <bits/stdc++.h>
using namespace std;
typedef double elm;
typedef long long lli;
typedef vector<elm> vec;
typedef vector<vec> mat;

const elm eps = 1e-8;

mat getE(int n) {
  mat e(n, vec(n, 0));
  for(int i = 0; i < n; ++i) {
    e[i][i] = 1;
  }
  return e;
}

mat multiply(const mat &a, const mat &b) {
  if(a.size() == 0 || b.size() == 0 ||
     a[0].size() != b.size()) return mat();
  int s = a.size();
  int t = b.size();
  int u = b[0].size();
  mat c(s, vec(u, 0));
  for(int i = 0; i < s; ++i) {
    for(int j = 0; j < u; ++j) {
      for(int k = 0; k < t; ++k) {
        c[i][j] += a[i][k] * b[k][j];
      }
    }
  }
  return c;
}

mat pow(const mat &a, const lli &n) {
  if(n == 0) return getE(a.size());
  mat b = pow(multiply(a, a), n/2);
  if(n & 1) return multiply(b, a);
  else return b;
}

// by Spaghetti Source
// O( n^3 ).
elm det(mat a) {
  const int n = a.size();
  elm res = 1;
  for (int i = 0; i < n; ++i) {
    int pivot = i;
    for (int j = i+1; j < n; ++j)
      if (abs(a[j][i]) > abs(a[pivot][i])) pivot = j;
    swap(a[pivot], a[i]);
    res *= a[i][i] * (i != pivot ? -1 : 1);
    if (abs(a[i][i]) < eps) break;
    for(int j = i+1; j < n; ++j)
      for(int k = n-1; k >= i; --k)
        a[j][k] -= a[i][k] * a[j][i] / a[i][i];
  }
  return res;
}

// by Spaghetti Source
// O( n^3 ).
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

int main() {

  return 0;
}
