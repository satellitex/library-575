#include<bits/stdc++.h>
using namespace std;

#define MAX_V 10000
int V;
vector<int> G[MAX_V];
int match[MAX_V];
bool used[MAX_V];


void graph_init(){
  for(int i=0;i<MAX_V;i++)G[i].clear();
}

void add_edge(int u,int v){
  G[u].push_back(v);
  G[v].push_back(u);
}
 
bool dfs(int v){
  used[v]=true;
  for(int i=0;i<(int)G[v].size();i++){
    int u=G[v][i],w=match[u];
    if( w<0 ||( !used[w] && dfs(w) )){
      match[v]=u;
      match[u]=v;
      return true;
    }
  }
  return false;
}
 
int bipartite_matching(){
  int res=0;
  memset(match,-1,sizeof(match));
  for(int v=0;v<V;v++){
    if(match[v]<0){
      memset(used,0,sizeof(used));
      if(dfs(v)){
        res++;
      }
    }
  }
  return res;
}
 
