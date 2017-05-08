#include <cmath>
#include <algorithm>
using namespace std;

const double INF = 1e+9;
const double EPS = 1e-8;

// 三次元
struct P3{
  double x,y,z;
  P3 operator + (const P3 &a) const{ return (P3){x+a.x, y+a.y, z+a.z}; }
  P3 operator - (const P3 &a) const{ return (P3){x-a.x, y-a.y, z-a.z}; }
  P3 operator * (const double &c) const{ return (P3){x*c, y*c, z*c}; }
};
P3 operator * (double c, const P3 &a){ return (P3){c * a.x, c * a.y, c * a.z}; }

struct L3{
  P3 a,b;
};

struct Sphere{ //球体
  P3 c;
  double r;

  double surface_area(){ return M_PI * r * r * 4.0;}
  double area(){ return M_PI * r * r * r * 4.0 / 3.0;}
};

double norm(P3 a){ return a.x*a.x + a.y*a.y + a.z*a.z; }

double abs(P3 a){ return sqrt(a.x*a.x + a.y*a.y + a.z*a.z); }

double dot(P3 a, P3 b){ return a.x*b.x+ a.y*b.y + a.z*b.z; }

P3 cross(P3 a, P3 b){
  return (P3){a.y*b.z - a.z*b.y, a.z*b.x - a.x*b.z, a.x*b.y - a.y*b.x};
}

bool cover(Sphere S, P3 p){
  return abs(S.c-p) < S.r + EPS;
}

double distPP(P3 a, P3 b){ return abs(a-b);};


// 線分と点の距離
double distSP(L3 S, P3 p){
  double a = dot(S.b-S.a, S.b-S.a);
  double b = dot(S.b-S.a, S.a- p );
  double c = dot(S.a- p , S.a- p );
  double t = -b / a;
  if(t < 0) return abs(S.a - p);
  if(t > 1) return abs(S.b - p);
  return a * t * t + 2 * b * t + c;
}

// 線分と線分の距離
double distSS(L3 S1, L3 S2){
  // a*s + b*t = e
  // c*s + d*t = f
  P3 v1 = S1.b-S1.a, v2 = S2.b-S2.a, v3 = S1.a-S2.a;
  double a = dot(v1,v1);
  double b = -dot(v1,v2);
  double c = -dot(v1,v2);
  double d = dot(v2,v2);
  double e = -dot(v1,v3);
  double f = dot(v2,v3);
  
  double ans = INF;
  ans = min(ans,distSP(S1,S2.a));
  ans = min(ans,distSP(S1,S2.b));
  ans = min(ans,distSP(S2,S1.a));
  ans = min(ans,distSP(S2,S1.b));

  double det=a*d-b*c;
  if(det != 0){
    double s = ( d*e-b*f) / det;
    double t = (-c*e+a*f) / det;
    if(0<=s && s<=1 && 0<=t && t<=1)
      ans = min(ans, distPP(S1.a+v1*s, S2.a+v2*t));
  }
  return ans;
}


/*
  極座標（半径,緯度,経度）から直交座標（x,y,z）を求める
  緯度・経度はラジアン
  r = sqrt(x*x + y*y + z*z)
  ido = acos(z / r) (0 <= ido <= pi)
  keido : cos(keido) = x / sqrt(x*x + y*y), sin(keido) = y / sqrt(x*x + y*y) を満たす実数 (0 <= keido < 2pi)
  未確認、wikipediaより引用
*/
P3 getP3(double r, double ido, double keido){
  P3 res;
  res.x = r * sin(ido) * cos(keido);
  res.y = r * sin(ido) * sin(keido);
  res.z = r * cos(ido);
  return res;
}


/*
  球面上にある二点間の球面距離
  緯度・経度はラジアン
*/
double gcdist(double r, double ido1, double keido1, double ido2, double keido2){
  double s = sin((ido1-ido2) / 2.0), t = sin((keido1-keido2) / 2.0);
  return 2.0 * r * asin(sqrt(max(0.0, s * s + cos(ido1) * cos(ido2) * t * t)));
}


/*
  球面上にある二点間の直線距離
  緯度・経度はラジアン
*/
double gcchord(double r, double ido1, double keido1, double ido2, double keido2){
  double c11 = cos(ido1), c12 = cos(keido1), c21 = cos(ido2), c22 = cos(keido2);
  double s11 = sin(ido1), s12 = sin(keido1), s21 = sin(ido2), s22 = sin(keido2);
  double dx = c11*c12 - c21*c22, dy = c11*s12 - c21*s22, dz = s11 - s21;
  return r * sqrt(dx*dx + dy*dy + dz*dz);
}


/* 
   4頂点A,B,C,Dからなる四面体の体積
   U : 辺AB の長さ
   V : 辺BC の長さ
   W : 辺CA の長さ
   u : 辺CD の長さ
   v : 辺AD の長さ
   w : 辺BD の長さ
   与えられた辺長から四面体が作れないときは, sqrt の引数が負になるかもしれない.
*/
double volume(double U, double V, double W, double u, double v, double w){
  double X = (w - U + v) * (U + v + w);
  double x = (U - v + w) * (v - w + U);
  double Y = (u - V + w) * (V + w + u);
  double y = (V - w + u) * (w - u + V);
  double Z = (v - W + u) * (W + u + v);
  double z = (W - u + v) * (u - v + W);
  double a = sqrt(x * Y * Z);
  double b = sqrt(y * Z * X);
  double c = sqrt(z * X * Y);
  double d = sqrt(x * y * z);
  return sqrt((-a + b + c + d) * (a - b + c + d) * (a + b - c + d) * (a + b + c - d))/(192 * u * v * w);
}


/*
  最小包含球
  二次元の場合は、Shere->Circle, P3->P にする。（元々二次元で書かれてあったが、三次元に拡張可能とのこと）
  Pは破壊される。
  スタックオーバーフローに注意。
  （未確認）
  平均計算量：O(n)
  最悪計算量：O(n^(DIM+1))
 */
class smallest_enclosing_ball{
  static const int DIM = 3; //三次元

  // 与えられた点集合
  int n;
  P3 *p;

  //境界にある点集合
  int n_bd;
  P3 bd[DIM+1];

  // {bd[0], ..., bd[n_bd-1]} の最小包含球を求める
  Sphere naive(){
    if(n_bd == 0) return (Sphere){(P3){0,0,0}, -1};

    // Gram-Schmidt で使う作業用
    P3 v[DIM]; // v[i] := ( v[0], ..., v[i-1] と直交するベクトル )
    double z[DIM]; // z[i] := |v[i]|ˆ2

    // (cen, sqrt(R2)) := ( {bd[0], ..., bd[i]} の最小包含球 )
    P3 cen = bd[0];
    double R2 = 0;
    for(int i=0;i<n_bd-1;i++){
      P3 delta = bd[i+1] - cen;

      // Gram-Schmidt
      v[i] = delta;
      for(int j=0;j<i;j++){
        v[i] = v[i] - dot(v[j], delta) / z[j] * v[j];
      }
      z[i] = dot(v[i], v[i]);

      double mag = norm(delta) - R2;
      cen = cen + mag / (2 * z[i]) * v[i];
      R2 += mag * mag / (4 * z[i]);
    }
    
    return (Sphere){cen, sqrt(R2)};
  }

  Sphere dfs(int i){
    if(i == n || n_bd == DIM + 1) return naive();
    
    // P[i] を境界として使わない場合
    Sphere c = dfs(i+1);
    // P[i] を境界として使う場合
    if(c.r==-1 || !cover(c,p[i])){
      bd[n_bd++] = p[i];
      c = dfs(i+1);
      n_bd--;
    }
    return c;
  }

public:
  Sphere solve(int n, P3 *p){
    this->n = n;
    this->p = p;
    random_shuffle(p, p+n);
    n_bd = 0;
    return dfs(0);
  }
};
