#include<bits/stdc++.h>
using namespace std;
typedef complex<double> P;
typedef complex<double> V;
typedef vector<P> vecP;
typedef pair<P,P> L;
typedef pair<P,P> S;
typedef pair<P,double> C;
const double eps=1e-8;
const double PI=acos(-1);
const double PI2=PI*2.0;
 
namespace std{
  bool operator < (const P &a,const P &b){
    return (a.real()!=b.real()?
            a.real() <b.real():
            a.imag() <b.imag());
  }
};
 
V normal(V a){
  assert( abs(a)>0 );
  return a/abs(a);
}
 
double Sqrt( double x ){
  if(x<0)return 0;
  else return sqrt(x);
}
 
P Vector(L a){
  return a.second-a.first;
}
 
bool eq(double a,double b){
  return (-eps<a-b&&a-b<eps);
}
 
bool eq(P a,P b){
  return ( eq(a.real(),b.real()) && eq(a.imag(),b.imag()) );
}
 
double dot(P a,P b){
  return real(b*conj(a));
}
 
double cross(P a,P b){
  return imag(b*conj(a));
}
 
double getArg(P a,P b){
  return arg(b*conj(a));
}
 
double getTime(V a,V b){
  assert( eq(cross(a,b),0) );
  return ( dot(a,b) < 0 ? -1.0 : 1.0 ) * abs(b) / abs(a);
}
 
 
P project(P a,P b,P c){
  b-=a,c-=a;
  return a+b*real(c/b);
}
 
P reflect(P a,P b,P c){
  b-=a,c-=a;
  return a+b*conj(c/b);
}
 
int ccw(P a,P b,P c){
  P ab=b-a,ac=c-a;
  P k=ac*conj(ab);
  if(k.imag()>eps)return 1;
  if(k.imag()<-eps)return -1;
  if(k.real()<-eps)return 2;
  if(abs(ab)+eps<abs(ac))return -2;
  return 0;
}
 
bool isParallel(P a,P b){
  return eq(0, cross(a,b));
}
 
bool isParallel(S a,S b){
  return eq(0, cross( Vector(a) , Vector(b) ) );
}
 
bool onLP(L l,P p){
  P a=l.first, b=l.second;
  return eq(0, cross(b-a,p-a));
}
 
bool onSP(S s,P p){
  P a=s.first, b=s.second;
  return eq( abs(b-a) , abs(a-p)+abs(b-p) );
}
 
bool isCrossSS(S s0,S s1){
  P a=s0.first, b=s0.second;
  P c=s1.first, d=s1.second;
  int f0 = ccw(a,b,c) * ccw(a,b,d);
  int f1 = ccw(c,d,a) * ccw(c,d,b);
  return (f0<=0 && f1<=0);
}
 
bool isCrossLS(L l,S s){
  P a=l.first, b=l.second;
  P c=s.first, d=s.second;
  return ( ccw(a,b,c) * ccw(a,b,d) <= 0 );
}
 
double distLP(L l,P p){
  P a=l.first, b=l.second;
  double res = cross(b-a,p-a) / abs(b-a);
  return abs(res);
}
 
double distSP(S s,P p){
  P a=s.first, b=s.second;
  if( dot(b-a,p-a) < eps )return abs(p-a);
  if( dot(a-b,p-b) < eps )return abs(p-b);
  return distLP(s,p);
}
 
double distSS(S s0,S s1){
  if( isCrossSS(s0,s1) )return 0;
  double res0 = min( distSP( s0, s1.first ) , distSP(s0, s1.second) );
  double res1 = min( distSP( s1, s0.first ) , distSP(s1, s0.second) );
  return min(res0,res1);
}
 
P getCrossLL(L l0,L l1){
  P a=l0.first, b=l0.second;
  P c=l1.first, d=l1.second;
  a-=d;b-=d;c-=d;
  return d+a+(b-a)*imag(a/c)/imag(a/c-b/c);
}
 
 
  
int inPolygon(vecP &t,P p){
  int n=t.size();
  double sum=0;
  for(int i=0;i<n;i++){
    P a=t[i],b=t[(i+1==n?0:i+1)];
    if( ccw(a,b,p)==0 )return 1;
    sum+= getArg(a-p,b-p);
  }
  if( abs(sum) < eps )return 0;
  else return 2;
}
 
vecP andrewScan(vecP &t){
  int N=t.size(),C=0;
  vecP R(N);
  for(int i=0;i<N;i++){
    while(2<=C&&ccw(R[C-2],R[C-1],t[i])==-1)C--;
    R[C++]=t[i];
  }
  vecP res(C);
  for(int i=0;i<C;i++)res[i]=R[i];
  return res;
}
  
vecP convexHull(vecP &t){
  sort(t.begin(),t.end());
  vecP u=andrewScan(t);
  reverse(t.begin(),t.end());
  vecP l=andrewScan(t);
  for(int i=1;i+1<(int)l.size();i++)u.push_back(l[i]);
  return u;
}
 
vecP cutConvex(vecP &t,L l){
  P a=l.first, b=l.second;
  int N=t.size();
  vecP res;
  for(int i=0;i<N;i++){
    P c=t[i],d=t[(i+1)%N];
    int C=ccw(a,b,c),D=ccw(a,b,d);
    if(C!=-1)res.push_back(c);
    if(C==-D&&abs(C)==1)res.push_back(getCrossLL( l ,L(c,d) ));
  }
  return res;
}
 
P getVector(const vecP &t, int id){
  int n=t.size();
  return t[ (id+1)%n ] - t[id%n];
}
 
double convex_diameter(vecP &t) {
  int n = t.size();
  int is = 0, js = 0;
  for (int i = 1; i < n; ++i) {
    if (imag(t[i]) > imag(t[is])) is = i;
    if (imag(t[i]) < imag(t[js])) js = i;
  }
  double maxd = norm(t[is]-t[js]);
  
  int i, maxi, j, maxj;
  i = maxi = is;
  j = maxj = js;
  do {
     
    if (cross( getVector(t,i), getVector(t,j)) >= 0) j = (j+1) % n;
     
    else i = (i+1) % n;
    if (norm(t[i]-t[j]) > maxd) {
      maxd = norm(t[i]-t[j]);
      maxi = i; maxj = j;
    }
  } while (i != is || j != js);
  return sqrt(maxd); /* farthest pair is (maxi, maxj). */
}
 
bool compare_y(const P &a,const P &b){
  return a.imag() < b.imag();
}
 
double closest_pair(P *a, int n){
  if(n <= 1) return 1e30;
  int m = n / 2;
  double x = a[m].real();
  double d = min(closest_pair(a, m), closest_pair(a + m, n - m));
  inplace_merge(a, a + m, a + n, compare_y);
  vector<P> b;
  for(int i=0;i<n;i++){
    if( abs(a[i].real() - x) >= d) continue;
    for(int j=0;j<(int)b.size();j++){
      double dx = real(a[i] - b[b.size() - j - 1]);
      double dy = imag(a[i] - b[b.size() - j - 1]);
      if(dy >= d) break;
      d = min(d, sqrt(dx * dx + dy * dy));
    }
    b.push_back(a[i]);
  }
  return d;
}
 
P _pool[200005];
double minDist(vecP &t){
  int n=t.size();
  for(int i=0;i<n;i++)_pool[i]=t[i];
  sort( _pool, _pool+n);
  return closest_pair(_pool, n);
}
 
int getStateCC(C a,C b){
  double ar=a.second, br=b.second;
  double dist=abs(a.first-b.first);
  if(dist>ar+br+eps)return 4;
  if(dist>ar+br-eps)return 3;
  if(dist>abs(ar-br)+eps)return 2;
  if(dist>abs(ar-br)-eps)return 1;
  return 0;
}
 
P getCrossCC(C a,C b){
  P p1=a.first, p2=b.first;
  double r1=a.second, r2=b.second;
  double cA = (r1*r1+norm(p1-p2)-r2*r2) / (2.0*r1*abs(p1-p2));
  return p1+(p2-p1)/abs(p1-p2)*r1*P(cA,Sqrt(1.0-cA*cA));
}
 
P getTangentCP_(C a,P p,int flg){
  P base=a.first-p;
  double ar=a.second;
  double w=Sqrt(norm(base)-ar*ar);
  P s=p+base*P(w,ar * flg)/norm(base)*w;
  return s;
}
 
vector<S> getTangentCP(C a,P p){
  vector<S> res;
  P s=getTangentCP_(a,p,1);
  P t=getTangentCP_(a,p,-1);
   
  if( eq(s,t) ){
    res.push_back( S(s, s+(a.first-p)*P(0,1) ) ); 
  }else{
    res.push_back( S(p,s) );
    res.push_back( S(p,t) );
  }
  return res;
}
 
S getInTangent(C a,C b,double flg=1.0){
  P ap=a.first,bp=b.first;
  double ar=a.second,br=b.second;
   
  P base=bp-ap;
  double w=ar+br;
  double h=Sqrt(norm(base)-w*w);
  P k=base*P(w,h*flg)/norm(base);
  return S(ap+k*ar,bp-k*br);
}
   
S getOutTangent(C a,C b,double flg=1.0){
  P ap=a.first,bp=b.first;
  double ar=a.second,br=b.second;
   
  P base=bp-ap;
  double h=br-ar;
   
  double w=Sqrt(norm(base)-h*h);
  P k=base*P(w,h*flg)/norm(base)*P(0,flg);
  return S(ap+k*ar,bp+k*br);
}
   
vector<S> getTangentCC(C a,C b){
  P ap=a.first,bp=b.first;
  double ar=a.second,br=b.second;
   
  vector<S> res;
  double dist=abs(ap-bp);
     
  if(dist>ar+br+eps)
    res.push_back(getInTangent(a,b,1));
   
  if(dist>ar+br-eps)
    res.push_back(getInTangent(a,b,-1));
   
  if(dist>abs(ar-br)+eps)
    res.push_back(getOutTangent(a,b,1));
   
  if(dist>abs(ar-br)-eps)
    res.push_back(getOutTangent(a,b,-1));
   
  return res;
}
 
 
vecP getCrossCL(C cir,L l){
  P a=l.first, b=l.second;
  double cr=cir.second;
  P cp=cir.first;
   
  vecP res;
  P base=b-a,  target=project(a,b,cp);
   
  double length=abs(base), h=abs(cp-target);
  base/=length;
   
  if(cr+eps<h)return res;
  double w=Sqrt(cr*cr-h*h);
  double L=getTime( normal(b-a) ,target-a)-w,  R=L+w*2.0;
   
  res.push_back(a+base*L);
  if( eq(L,R) )return res;
  res.push_back(a+base*R);
 
  return res;
}
 
vecP getCrossCS(C cir,S s){
  vecP tmp=getCrossCL(cir,s);
  vecP res;
  for(int i=0;i<(int)tmp.size();i++)
    //  if( ccw(s.first,s.second, tmp[i] ) == 0)
    if( onSP( s, tmp[i] ) )
      res.push_back(tmp[i]);
  return res;
}
 
double getArea(C c,P a,P b){
  P cp=c.first;
  double cr=c.second;
   
  P va=cp-a,  vb=cp-b;
  double A=abs(va), B=abs(vb);
  double f=cross(va,vb), d=distSP( S(a,b) ,cp), res=0;
   
  if( eq(0, f ) )return 0;
  if(A<cr+eps&&B<cr+eps)return f*0.5;
  if(d>cr-eps)return cr*cr*PI*getArg(va,vb)/PI2;
    
  vecP u=getCrossCS(c, S(a,b) );
   
  assert( !u.empty() );
  u.insert(u.begin(), a),  u.push_back(b);
  
  for(int i=0;i+1<(int)u.size();i++) res+=getArea(c,u[i],u[i+1]);
  return res;
}
  
double getCrossArea(vecP t,C c){
  int n=t.size();
  if(n<3)return 0;
  double res=0;
  for(int i=0;i<n;i++){
    P a=t[i], b=t[(i+1)%n];
    res+=getArea(c,a,b);
  }
  return res;
}
 
 
double calcPolygonArea(const vecP &t){
  double res=0;
  int n=t.size();
  for(int i=0;i<n;i++){
    res+= cross( t[ (i+1)%n ],t[i] );
  }
  return abs(res)*0.5;
}
 
P inputP(){
  int x,y;
  scanf("%d %d",&x,&y);
  return P(x,y);
}
 
 
void pr(P p,string str="\n"){
  printf("%.10f %.10f",p.real(),p.imag());
  cout<<str;
}
 
 
struct edge{
  int to;
  double cost;
};
 
typedef vector<edge> Node;
typedef vector<Node> Graph;
 
bool dfs(int pos,double sum,Graph &G, vector<int> &visited , vector< double > &mem){
  if(visited[pos]){
    return !eq( sum-mem[pos] , 0 );
  }
   
  visited[pos]=true;
  mem[pos]=sum;
   
  for(int i=0;i<(int)G[pos].size();i++){
    edge e=G[pos][i];
    if( dfs(e.to, sum+e.cost, G, visited, mem) )return true;
  }
  return false;
}
 
bool solve(Graph &G){
  int N=G.size();
  vector<int> visited(N,false);
  vector<double> mem(N);
   
  for(int i=0;i<N;i++){
    if(visited[i])continue;
    if(dfs(i,0,G,visited,mem))return true;
  }
   
  return false;
}
 
Graph segments2graph( vector<S> segments ){
  Graph G;
  vector< P > points;
  int n=segments.size();
  for(int i=0;i<n;i++){
    S si=segments[i];
    points.push_back( si.first );
    points.push_back( si.second );
    for(int j=0;j<i;j++){
      S sj=segments[j];
      if( isCrossSS( si , sj ) && !isParallel( si, sj ) ){
        points.push_back( getCrossLL(si,sj) );
      }
    }
  }
 
  sort(points.begin(), points.end());
  points.erase( unique( points.begin(), points.end() ) , points.end() );
   
  G.resize(points.size());
 
  typedef pair< double , int > Pair;
     
  for(int i=0;i<n;i++){
    S si=segments[i];
    P sf=si.first, se=si.second;
    vector< Pair > targets;
    for(int j=0;j<(int)points.size();j++){
      P pj=points[j];
      // if( ccw( si.first, si.second, pj ) == 0 ){
      if( onSP( si, pj ) ){
        targets.push_back( Pair(getTime(se-sf ,pj-sf ), j ) );
      }
    }
    /*
    for(int j=0;j<(int)targets.size();j++){
      for(int k=0;k<j;k++){
        Pair a=targets[j];
        Pair b=targets[k];
        P ap=points[ a.second ];
        P bp=points[ b.second ];
        G[a.second].push_back( (edge){ b.second , getArg(ap,bp) } );
        G[b.second].push_back( (edge){ a.second , getArg(bp,ap) } );
      }
    }
    */
 
    sort( targets.begin(), targets.end() );
    for(int j=0;j+1<(int)targets.size();j++){
      Pair a=targets[j];
      Pair b=targets[j+1];
      P ap=points[ a.second ];
      P bp=points[ b.second ];
      G[a.second].push_back( (edge){ b.second , getArg(ap,bp) } );
      G[b.second].push_back( (edge){ a.second , getArg(bp,ap) } );
    }
 
  }
 
  return G;
}
 
int main(){
  while(1){
    int n;
    scanf("%d",&n);
    if(n==0)break;
    vector< S > segments;
    for(int i=0;i<n;i++){
      P a=inputP();
      P b=inputP();
      segments.push_back( S(a,b) );
    }
    Graph G=segments2graph(segments);
 
 
    if( solve(G) )cout<<"yes"<<endl;
    else cout<<"no"<<endl;
  }
  return 0;
}
