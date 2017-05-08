#include <bits/stdc++.h>
using namespace std;

typedef complex<double> P;

const double eps = 1e-8;
const double inf = 1e20;

double equals(double a, double b) {
  return abs(a-b) < eps;
}

double dot(P a, P b) {
  return a.real()*b.real() + a.imag()*b.imag();
}

// ax^2 + bx + c = 0
vector<double> quadraticFormula(double a, double b, double c) {
  vector<double> x;
  if(a == 0.0) {
    if(b != 0.0) x.push_back(-c/b);
  } else {
    double d = b*b - 4.0*a*c;
    if(d == 0.0) {
      x.push_back(-b/(2*a));
    } else if(d > 0.0) {
      x.push_back((-b-sqrt(d))/(2*a));
      x.push_back((-b+sqrt(d))/(2*a));
    }
  }
  return x;
}

// ax^2 + 2hbx + c = 0
vector<double> quadraticFormula2(double a, double hb, double c) {
  vector<double> x;
  if(a == 0.0) {
    if(hb != 0.0) x.push_back(-c/(2.0*hb));
  } else {
    double d = hb*hb - a*c;
    if(d == 0.0) {
      x.push_back(-hb/a);
    } else if(d > 0.0) {
      x.push_back((-hb-sqrt(d))/a);
      x.push_back((-hb+sqrt(d))/a);
    }
  }
  return x;
}

/*
  2つの円（超球）が等速直線運動し、半径が線形的に変化するとき、
  2つの円が衝突する時間を返す。Pはn次元に対応。
 */
double calcT(P p1, P v1, double r1, double vr1,
             P p2, P v2, double r2, double vr2) {
  double a = norm(v1 - v2) - pow(vr1 + vr2, 2);
  double hb = dot(v1, p1) - dot(v1, p2) - dot(v2, p1) + dot(v2, p2) - (r1 + r2)*(vr1 + vr2);
  double c = norm(p1 - p2) - pow(r1 + r2, 2);
  vector<double> v = quadraticFormula2(a, hb, c);
  for(int i = 0; i < v.size(); ++i) {
    if(v[i] >= 0) return v[i];
  }
  return inf;
}

bool hit(P p1, double r1, P p2, double r2) {
  double d = abs(p1 - p2);
  return equals(d, r1 + r2) || d < r1 + r2;
}

// AOJ 1139
int main() {
  int N;
  double T, R;
  for(; cin >> N >> T >> R && (N|(T!=0)|(R!=0));) {
    R /= 2.0;
    vector<string> name(N);
    vector<P> ps(N);
    vector<deque<double> > ts(N);
    vector<deque<P> > vs(N);
    for(int i = 0; i < N; ++i) {
      double t;
      P v;
      cin >> name[i];
      cin >> t >> ps[i].real() >> ps[i].imag();
      while(1) {
        cin >> t >> v.real() >> v.imag();
        ts[i].push_back(t);
        vs[i].push_back(v);
        if(t == T) break;
      }
    }

    vector<int> flag(N, false);
    flag[0] = true;
    for(int k = 0; k < N; ++k) {
      for(int i = 0; i < N; ++i) {
        for(int j = i+1; j < N; ++j) {
          if(hit(ps[i], R, ps[j], R) && (flag[i] || flag[j])) {
            flag[i] = flag[j] = true;
          }
        }
      }
    }

    double now = 0;
    while(!equals(now, T)) {
      double dt = inf;
      for(int i = 0; i < N; ++i) {
        dt = min(dt, ts[i].front() - now);
      }
      for(int i = 0; i < N; ++i) {
        if(!flag[i]) continue;
        for(int j = 0; j < N; ++j) {
          if(flag[j]) continue;
          dt = min(dt, calcT(ps[i], vs[i].front(), R, 0, ps[j], vs[j].front(), R, 0));
        }
      }
      if(dt == inf) break;
      now += dt;
      for(int i = 0; i < N; ++i) {
        ps[i] += dt * vs[i].front();
        if(equals(now, ts[i].front())) {
          vs[i].pop_front();
          ts[i].pop_front();
        }
      }

      for(int i = 0; i < N; ++i) {
        for(int j = i+1; j < N; ++j) {
          if(hit(ps[i], R, ps[j], R) && (flag[i] || flag[j])) {
            flag[i] = flag[j] = true;
          }
        }
      }
    }

    vector<string> res;
    for(int i = 0; i < N; ++i) {
      if(flag[i]) res.push_back(name[i]);
    }
    sort(res.begin(), res.end());
    for(int i = 0; i < res.size(); ++i) {
      cout << res[i] << endl;
    }
  }
  return 0;
}
