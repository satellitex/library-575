#include<bits/stdc++.h>
using namespace std;
int INF=1e9;
 
struct UnionFind{
  vector<int> par,rank;
  void init(int n){
    par.clear();
    rank.clear();
    par.resize(n);
    rank.resize(n);
    for(int i=0;i<n;i++){
      par[i]=i;
      rank[i]=1;
    }
  }
   
  int find(int x){
    if(x==par[x])return x;
    return par[x]=find(par[x]);
  }
 
  bool same(int x,int y){
    return ( find(x)==find(y) );
  }
   
  void unite(int x,int y){
    x=find(x);
    y=find(y);
    if(x==y)return;
    if(rank[x]<rank[y])swap(x,y);
    par[y]=x;
    rank[x]+=rank[y];
  }
};
 
struct edge{
  int from,to,cost,id;
  bool operator < (const edge e)const{
    return cost > e.cost;
  }
};
 
typedef priority_queue< edge > prque;
typedef prque* Prque;
 
typedef vector<edge> Node;
typedef vector<Node> Graph;
 
 
 
 
 
Prque Merge(vector< Prque > &Q, vector<edge> &ev,int A,int C){
  if( Q[C]->size() < Q[A]->size() ){
 
    while( !Q[C]->empty() ){
      edge e=Q[C]->top();
      e.cost-=ev[C].cost;
      e.cost+=ev[A].cost;
      //      cout<<e.from<<' '<<e.to<<' '<<e.cost<<' '<<ev[A].cost<<' '<<ev[C].cost<<endl;
      Q[A]->push(e);
      Q[C]->pop();
    }
    ev[C].cost=ev[A].cost;
    return Q[A];
  }else{
    while( !Q[A]->empty() ){
      edge e=Q[A]->top();
      e.cost-=ev[A].cost;
      e.cost+=ev[C].cost;
      //      cout<<e.from<<' '<<e.to<<' '<<e.cost<<endl;
      Q[C]->push(e);
      Q[A]->pop();
    }
    return Q[C];
  }
}
 
int solve(Graph &G,vector<edge> &edges,int root){
 
  int n=G.size(), res=0;
  vector<int> used(n,0);
  vector< edge > ev(n, (edge){0,0,0,-1} );
  vector< prque > pool(n);
  vector< Prque > Q(n);
  for(int i=0;i<n;i++)Q[i]=&pool[i];
   
  UnionFind uf;
  uf.init(n);
   
  for(int i=0;i<(int)edges.size();i++){
    edge e=edges[i];
    Q[ e.to ]->push( e );
  }
   
  used[root]=2;
  for(int Pos=0;Pos<n;Pos++){
    if(used[Pos]==2)continue;
    int pos=Pos;
    vector<int> path;
     
    while( used[pos] != 2 ){
      pos=uf.find(pos);
       
      used[pos]=1;
      path.push_back(pos);
      if( Q[pos]->empty() ){
        return INF;
      }
       
      edge e=Q[pos]->top();
 
       
      Q[pos]->pop();
      e.cost-=ev[pos].cost;
      if( uf.same(e.from,pos) ) continue;
      int tmpcost=ev[pos].cost;
      /*
      cout<<" pos="<<pos;
      cout<<" e.from="<<e.from;
      cout<<" e.to="<<e.to;
      cout<<" e.cost="<<e.cost;
      cout<<" tmpcost="<<tmpcost<<endl;
      cout<<endl;
      */
      res+=e.cost;
      e.cost+=tmpcost;
      ev[pos]=e;
      if( used[ uf.find(e.from) ] == 2 )break;
      if( used[ uf.find(e.from) ] == 0 ){
        pos=e.from;
        continue;
      }
      int pre=uf.find(e.from);
      set<int> mp;
      while(1){
        if(mp.count(pre))break;
        mp.insert(pre);
        
        if(!uf.same(pre,pos)){
          int A=uf.find(pre), B=uf.find(pos);
          uf.unite(A,B);
          int C=uf.find(A);
          /*
          cout<<" !!A="<<A;
          cout<<" !!B="<<B;
          cout<<" !!C="<<C<<endl;
          */
          Prque tmp=NULL;
          if(B==C)tmp=Merge(Q,ev,A,C);
          else if(A==C)tmp=Merge(Q,ev,B,C);
          else assert(0);
           
          Q[C]=tmp;
        }
        pre=uf.find(ev[pre].from);
      }
    }// while_pos
 
    for(int i=0;i<(int)path.size();i++)used[ path[i] ]=2;
  }// Pos
  return res;
}
 
int main(){
  int V,E,r;
  vector<edge> edges;
  Graph G;
  cin>>V>>E>>r;
  G.resize(V);
  for(int i=0;i<E;i++){
    int a,b,c;
    cin>>a>>b>>c;
    G[a].push_back( (edge){a,b,c,i} );
    edges.push_back( (edge){a,b,c,i} );
  }
  cout<< solve( G, edges, r ) <<endl;  
  return 0;
}
