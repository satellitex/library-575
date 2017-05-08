// last update 2016/05/03

#include <bits/stdc++.h>
using namespace std;

typedef complex<double> P;
typedef vector<P> G;

const double eps = 1e-8;

bool equals(double a, double b) { return abs(a - b) < eps; }

bool equals(P a, P b) {
  return equals(a.real(), b.real()) && equals(a.imag(), b.real());
}

namespace std {
  bool operator < (const P &a, const P b) {
    return ( a.real() != b.real()
             ? a.real() < b.real() : a.imag() < b.imag() );
  }
}

struct L : public vector<P> {
  L(P a, P b) {
    push_back(a), push_back(b);
  }
};

struct C {
  P p;
  double r;
  C() {}
  C(P p, double r) : p(p), r(r) {}
};

double dot  (P a, P b) { return a.real()*b.real() + a.imag()*b.imag(); }
double cross(P a, P b) { return a.real()*b.imag() - a.imag()*b.real(); }

bool isOrthogonal(P a, P b) { return equals(dot(a, b), 0.0); }
bool isParallel  (P a, P b) { return equals(cross(a, b), 0.0); }

P project(P s1, P s2, P p) {
  P base = s2 - s1;
  double t = dot(p - s1, base)/norm(base);
  return s1 + base*t;
}

P reflect(P s1, P s2, P p) { return p + (project(s1, s2, p) - p)*2.0; }

P getRotateP(P p, double rad, P o = P(0, 0)) {
  P q = p - o;
  return o + P(q.real()*cos(rad) - q.imag()*sin(rad),
               q.real()*sin(rad) + q.imag()*cos(rad));
}

namespace CCW {
  enum { COUNTER_CLOCKWISE = 1, CLOCKWISE = -1,
         ONLINE_BACK = 2, ONLINE_FRONT = -2,
         ONSEGMENT = 0 };
}

int ccw(P p0, P a, P b) {
  a -= p0;
  b -= p0;
  if(cross(a, b) > eps) return 1;
  if(cross(a, b) < -eps) return -1;
  if(dot(a, b) < -eps) return 2;
  if(norm(b)-norm(a) > eps) return -2;
  return 0;
}

// 0:外部, 1:境界, 2:内部
int contain(C c, P p) {
  if(abs(c.p - p) - c.r <  -eps) return 2;
  if(abs(c.p - p) - c.r >  eps) return 0;
  return 1;
}

// 点が多角形の内部/境界/外部のどこにあるかを判定する．
// 0:外部, 1:境界, 2:内部
// by Spaghetti Source
int contain(const vector<P> &g, const P &p) {
  int n = g.size();
  int in = 0;
  for(int i = 0; i < n; ++i) {
    P a = g[i] - p, b = g[(i+1)%n] - p;
    if(a.imag() > b.imag()) swap(a, b);
    if(a.imag() < eps && eps < b.imag() && cross(a, b) < -eps) in = !in;
    if(abs(cross(a, b)) < eps && dot(a, b) < eps) return 1;
  }
  return in * 2;
}

// 線分と線分の交差判定。
// T字やL字みたいになっている時もtrueを返すが、
// これを除外したい場合に <= を < に変更してもうまくいかない。
bool isIntersect(P a1, P a2, P b1, P b2) {
  return ( ccw(a1, a2, b1) * ccw(a1, a2, b2) <= 0 &&
           ccw(b1, b2, a1) * ccw(b1, b2, a2) <= 0 );
}

double getDistanceLP(P s1, P s2, P p) {
  return abs(cross(s2 - s1, p - s1)/abs(s2 - s1));
}

double getDistanceSP(P s1, P s2, P p) {
  if(dot(s2 - s1, p - s1) < 0.0) return abs(p - s1);
  if(dot(s1 - s2, p - s2) < 0.0) return abs(p - s2);
  return getDistanceLP(s1, s2, p);
}

double getDistance(P a1, P a2, P b1, P b2) {
  if(isIntersect(a1,a2,b1,b2)) return 0.0;
  return min(min(getDistanceSP(a1,a2,b1), getDistanceSP(a1,a2,b2)),
             min(getDistanceSP(b1,b2,a1), getDistanceSP(b1,b2,a2)));
}

// 直線と円の交差判定。戻り値は交点の数。
int isIntersect(P s1, P s2, C c) {
  double d = getDistanceLP(s1, s2, c.p);
  if(equals(d, c.r)) return 1;
  else if(d < c.r) return 2;
  else return 0;
}

// 円と円の交差判定。
//  0 : 交差、内包なし
//  1 : 外部で1点と接する
//  2 : 2点で交差
// -1 : 内包して接する
// -2 : 完全に内包
int isIntersect(C a, C b) {
  double x = a.p.real() - b.p.real();
  double y = a.p.imag() - b.p.imag();
  double s = a.r + b.r;
  double d = x*x + y*y;
  s *= s;
  if(equals(d, s)) return 1;
  if(d > s) return 0;
  double r = abs(a.r - b.r);
  r *= r;
  if(equals(d, r)) return -1;
  if(d > r) return 2;
  return -2;
}

// 直線と直線の交点。
P getCrossP(P a1, P a2, P b1, P b2) {
  P a = a2 - a1;
  P b = b2 - b1;
  // cross の符号関係あり
  return a1 + a * cross(b, b1 - a1)/cross(b, a);
}

// 直線と円の交点。
vector<P> getCrossPLC(P s1, P s2, C c) {
  vector<P> v;
  P p = project(s1, s2, c.p);
  double dist = getDistanceLP(s1, s2, c.p);
  if (equals(dist, c.r)) {
    v.push_back(p);
  } else if (dist < c.r) {
    double h = abs(p - c.p);
    double d = sqrt(c.r * c.r - h * h);
    P base = s2 - s1;
    // push_backする順番は重要
    v.push_back(p - d * base / abs(base));
    v.push_back(p + d * base / abs(base));
  }
  return v;
}

// 円と円の交点。
vector<P> getCrossP(C c1, C c2) {
  vector<P> v;
  int cp = isIntersect(c1,c2);
  if(cp == 0 || cp == -2) return v;
  
  double ll = norm(c1.p - c2.p);
  double A = ( c1.r * c1.r - c2.r * c2.r + ll ) / ( 2.0 * ll );
  P base = c2.p - c1.p;

  if(abs(cp) == 1) {
    v.push_back(c1.p + A*base);
  } else {
    P n(-base.imag(), base.real());
    n /= abs(n);
    double h = sqrt(c1.r * c1.r - A*A*ll);
    v.push_back(c1.p + A*base + h*n);
    v.push_back(c1.p + A*base - h*n);
  }
  return v;
}

// 2つのベクトルのなす角(0 <= rad <= PI)
// AOJ2233でバグったので更新(2015/05/29)
double getAngle(P a, P b) {
  double v = dot(a, b) / (abs(a) * abs(b));
  if(v > 1.0) return 0;
  if(v < -1.0) return M_PI;
  return acos(v);
}

// a -> b への角度の移動量(0 <= rad < 2*PI)
double getAngleVector(P a, P b) {
  double A = arg(a);
  double B = arg(b);
  double rad = B-A;
  rad = fmod(rad, M_PI*2.0);
  if(rad < 0.0) rad += M_PI*2.0;
  return rad;
}

// 2つの直線に接する半径rの円
C getC(P a1, P a2, P b1, P b2, double r) {
  P a = a2 - a1;
  P b = b2 - b1;
  if(cross(a, b) < 0) swap(a, b);
  P p = getCrossP(a1, a2, b1, b2);
  double rad = getAngle(a, b);
  double alpha = arg(a);
  double d = r/sin(rad/2.0);
  C res;
  res.r = r;
  res.p = p + d * P(cos(alpha + rad/2.0), sin(alpha + rad/2.0));
  return res;
}

// ヘロンの公式を用いて三角形の3辺の長さから面積を求める
double heron(double a, double b, double c) {
  double s = (a+b+c)/2.0;
  return sqrt(s*(s-a)*(s-b)*(s-c));
}

// 多角形の面積(符号付き)
double getArea(vector<P> &G) {
  int n = G.size();
  double S = 0;
  for(int i = 0; i < n; ++i) {
    S += cross(G[i], G[(i+1)%n]);
  }
  return S/2.0;
}

double getTriArea(P a, P b, P c) {
  return cross(b-a, c-a) / 2.0;
}

double getSectorArea(P a, P b, C c) {
  return getAngleVector(a-c.p, b-c.p) * c.r * c.r / 2.0;
}

double getBowArea(P a, P b, C c) { // a,b is a chord of c.
  return getSectorArea(a, b, c) - getTriArea(a, b, c.p);
}

double getArea(P a, P b, C c) {
  int d = 1;
  if(cross(a,b) < 0.0) {
    d = -1;
    swap(a, b);
  }
  if(cross(a, b) < eps) return 0;

  G t(3);
  t[0] = P(0, 0), t[1] = a, t[2] = b;

  if(contain(c, t[0]) && contain(c, t[1]) && contain(c, t[2])) {
    return d * getTriArea(t[0], t[1], t[2]);
  }

  if(getDistanceSP(t[0], t[1], c.p) - c.r > -eps &&
     getDistanceSP(t[1], t[2], c.p) - c.r > -eps &&
     getDistanceSP(t[2], t[0], c.p) - c.r > -eps) {
    if(cross(t[1] - t[0], c.p - t[0]) > -eps &&
       cross(t[2] - t[1], c.p - t[1]) > -eps &&
       cross(t[0] - t[2], c.p - t[2]) > -eps) {
      return d * (c.r * c.r * M_PI);
    } else {
      return 0;
    }
  }

  double S = c.r * c.r * M_PI;
  vector<P> cp(6);
  vector<int> exist(6);
  for(int i = 0; i < 3; ++i) {
    P p = t[i], q = t[(i+1)%3];
    if(getDistanceLP(p, q, c.p) - c.r > -eps) {
      if(cross(q - p, c.p - p) < 0) return 0;
      continue;
    }
    if(getDistanceSP(p, q, c.p) - c.r > -eps) continue;
    vector<P> v = getCrossPLC(t[i], t[(i+1)%3], c);
    if(v.size() == 2 && !equals(v[0], v[1])) {
      S -= getBowArea(v[0], v[1], c);
      for(int j = 0; j < v.size(); ++j) {
        int k = (i*2 + j*3) % 6;
        cp[k] = v[j];
        exist[k] = true;
      }
    }
  }
  for(int i = 0; i < 3; ++i) {
    int a = i*2, b = i*2+1;
    if(!exist[a] || !exist[b]) continue;
    if(contain(c, t[i]) == 2 && !equals(cp[a], cp[b])) {
      S += getBowArea(cp[a], cp[b], c) + getTriArea(t[i], cp[a], cp[b]);
    }
  }
  return S * d;
}

// 多角形と円の共通部分の面積(符号付き)
double getArea(G g, C c) {
  int n = g.size();
  double res = 0;
  for(int i = 0; i < n; ++i) {
    res += getArea(g[i], g[(i+1)%n], c);
  }
  return res;
}

// 三角形の内接円の半径
double getIncircleR(P p1, P p2, P p3) {
  double a = abs(p1 - p2);
  double b = abs(p2 - p3);
  double c = abs(p3 - p1);
  return heron(a,b,c)*2.0/(a+b+c);
}

vector<P> convex_hull(vector<P> ps) {
  int n = ps.size(), k = 0;
  sort(ps.begin(), ps.end());
  vector<P> ch(2*n);
  for(int i = 0; i < n; ch[k++] = ps[i++])
    while(k >= 2 && ccw(ch[k-2], ch[k-1], ps[i]) <= 0) --k;
  for(int i = n-2, t = k+1; i >= 0; ch[k++] = ps[i--])
    while(k >= t && ccw(ch[k-2], ch[k-1], ps[i]) <= 0) --k;
  ch.resize(k-1);
  return ch;
}

vector<P> convex_cut(vector<P> ps, P p1, P p2, int dir = CCW::CLOCKWISE) {
  vector<P> v;
  for(int i = 0; i < ps.size(); ++i) {
    P a = ps[i];
    P b = ps[(i+1)%ps.size()];
    if(ccw(p1, p2, a) != dir) v.push_back(a);
    if(ccw(p1, p2, a)*ccw(p1, p2, b) == -1)
      v.push_back(getCrossP(a, b, p1, p2));
  }
  return v;
}

// キャリパー法を用いて多角形の最も遠い２頂点の距離を返す
// O(NlogN)
// by Ari book
double caliper(const vector<P> &ps) {
  vector<P> qs = convex_hull(ps);
  int n = qs.size();
  if(n == 2) {
    return abs(qs[0] - qs[1]);
  }
  int i = 0, j = 0;
  for(int k = 0; k < n; ++k) {
    if(!(qs[i] < qs[k])) i = k;
    if(qs[j] < qs[k]) j = k;
  }
  double res = 0;
  int si = i, sj = j;
  while(i != sj || j != si) {
    res = max(res, abs(qs[i] - qs[j]));
    if(cross((qs[(i+1)%n] - qs[i]), (qs[(j+1)%n] - qs[j])) < 0) {
      i = (i+1)%n;
    } else {
      j = (j+1)%n;
    }
  }
  return res;
}

// 円と円の共通接線
// validate by AOJ2201 2016/05/03, RUPC2016day2I 2016/05/03
// 戻り値は外接線2本、内接線2本の順で並び、
// 各直線の向きは小さい円->大きい円
// それぞれ、2本がほぼ同じ直線であれば1本に集約され、存在しなければ0本になる。
vector<L> getLine(C a, C b) {
  vector<L> res;
  if (a.r > b.r) swap(a, b);
  if (a.p == b.p) return res;
  P base = b.p - a.p;
  for (double s = -1; s <= 1; ++++s) {
    double c, d, ee, e;
    ee = norm(base);
    e = sqrt(ee);
    d = b.r + a.r * s;
    if (abs(e - d) < eps) {
      P v = base * P(0, 1);
      P m = b.p - base * b.r / abs(base);
      res.push_back(L(m, m + v));
    } else if (ee - d * d >= 0) {
      c = sqrt(ee - d * d);
      for (double t = -1; t <= 1; ++++t) {
        P rotate(c / e, d / e * t);
        P v = -base * rotate * s;
        P nv = v / abs(v) * P(0, t);
        res.push_back(L(a.p + nv * a.r, b.p - nv * b.r * s));
      }
    }
  }
  return res;
}

// 凸多角形の重心を求める。
// verified by UVA 10002 2013/07/07
P getCenterP(const vector<P> &g) {
  int n = g.size();
  P res(0,0);
  double sum = 0;
  for(int i = 1; i+1 < n; ++i) {
    double s = heron(abs(g[i]-g[0]),abs(g[i+1]-g[i]),abs(g[0]-g[i+1]));
    res += (g[0]+g[i]+g[i+1])*s/3.0;
    sum += s;
  }
  return res/sum;
}

// verified by AOJ 2454 2015/05/14
struct Edge { int to; double w; };
typedef vector<vector<Edge> > Graph;
Graph segmentArrangement(const vector<P> &s1, const vector<P> &s2,
                         vector<P> &ps) {
  int m = s1.size();
  for(int i = 0; i < m; ++i) {
    ps.push_back(s1[i]);
    ps.push_back(s2[i]);
    for(int j = i+1; j < m; ++j) {
      if(equals(cross(s1[i] - s2[i], s1[j] - s2[j]), 0.0)) continue;
      if(!isIntersect(s1[i], s2[i], s1[j], s2[j])) continue;
      ps.push_back(getCrossP(s1[i], s2[i], s1[j], s2[j]));
    }
  }
  sort(ps.begin(), ps.end());
  ps.erase(unique(ps.begin(), ps.end()), ps.end()); // eps比較しなくてええの？

  int n = ps.size();
  Graph g(n);
  for(int i = 0; i < m; ++i) {
    vector<pair<double, int> > v;
    for(int j = 0; j < n; ++j) {
      if(!ccw(s1[i], s2[i], ps[j])) {
        v.push_back(make_pair(norm(s1[i] - ps[j]), j));
      }
    }
    sort(v.begin(), v.end());
    for(int j = 0; j+1 < v.size(); ++j) {
      int a = v[j].second, b = v[j+1].second;
      double w = abs(ps[a] - ps[b]);
      g[a].push_back((Edge){b, w});
      g[b].push_back((Edge){a, w});
    }
  }
  return g;
}

int main(void) {
}

