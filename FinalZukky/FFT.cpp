/*
  FFT from Spaghetti Source
 */
#include <bits/stdc++.h>
using namespace std;
typedef complex<double> P;

const double pi = acos(-1);

// O(NlogN)
// from Spaghetti Source
vector<P> fft(vector<P> a, bool inv = false) {
  int n = a.size();
  double theta = 2.0 * pi / n;
  if (inv) theta *= -1;

  for (int m = n; m >= 2; m >>= 1) {
    int mh = m >> 1;
    for (int i = 0; i < mh; i++) {
      P w = exp(i*theta*P(0,1));
      for (int j = i; j < n; j += m) {
        int k = j + mh;
        P x = a[j] - a[k];
        a[j] += a[k];
        a[k] = w * x;
      }
    }
    theta *= 2;
  }
  int i = 0;
  for (int j = 1; j < n - 1; j++) {
    for (int k = n >> 1; k > (i ^= k); k >>= 1);
    if (j < i) swap(a[i], a[j]);
  }

  if (inv) for (int i = 0; i < n; ++i) a[i] /= n;

  return a;
}

/*
  G(x) = Σg[i]*x^i, H(x) = Σh[j]*x^j のとき、
  (G * H)(x) = G(x)*H(x) = Σgh[k]*x^k におけるgh[k]を返す。
*/
vector<P> multiply(vector<P> g, vector<P> h) {
  int n = 1;
  while (!(n > (int)g.size() - 1 + (int)h.size() - 1)) n *= 2;
  g.resize(n);
  h.resize(n);
  vector<P> gg = fft(g);
  vector<P> hh = fft(h);
  vector<P> ff(n);
  for (int i = 0; i < n; ++i) {
    ff[i] = gg[i] * hh[i];
  }
  vector<P> gh = fft(ff, true);
  return gh;
}

/*
  A[i] : 価値iの品物がA[i]種類
  B[j] : 価値jの品物がB[j]種類
  A,Bそれぞれから品物を１種類選び、
  価値の和がちょうどkとなるような場合の数C[k]を返す。
 */
vector<double> convolution(const vector<int> &A, const vector<int> &B) {
  int N = A.size() - 1;
  vector<P> g(N+1), h(N+1);
  for (int i = 0; i <= N; ++i) {
    g[i].real() = A[i];
    h[i].real() = B[i];
  }
  vector<P> gh = multiply(g, h);
  vector<double> C(2*N+1);
  for (int i = 0; i <= 2*N; ++i) {
    C[i] = round(gh[i].real());
  }
  return C;
}

/*
  AtCoder Typical Contest 001  C - 高速フーリエ変換

  主菜は、価格が i 円のものが Ai 種類あります (1≦i≦N)。
  副菜は、価格が i 円のものが Bi 種類あります (1≦i≦N)。
  定食は、主菜と副菜を 1 種類ずつ選んで構成します。 定食の価格は、選んだ主菜と副菜の価格の和とします。
  
  各 k(1≦k≦2N) について、価格が k 円になる定食の種類の数を計算して下さい。
*/
int main() {
  int N;
  cin >> N;
  vector<int> A(N+1), B(N+1);
  for(int i = 1; i <= N; ++i) {
    cin >> A[i] >> B[i];
  }
  vector<double> C = convolution(A, B);
  for(int i = 1; i <= 2*N; ++i) {
    cout << (int)C[i] << endl;
  }
  return 0;
}
