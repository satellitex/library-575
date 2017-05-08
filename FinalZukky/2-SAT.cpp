typedef pair<int,int> P;

//強連結成分分解こっから-----------------------------------
vector<int> G[20002];
vector<int> rG[20002];
bool used[20002];
int cmp[20002];
vector<int> vs;

void init(int size){
  vs.clear();
  for(int i=0;i<size;i++) {
    G[i].clear();
    rG[i].clear();
    cmp[i] = -1;
    used[i] = 0;
  }
}

void add_edge(int a,int b){//a->b
  G[a].push_back( b );
  rG[b].push_back( a );
}
void dfs(int id){
  if( used[id] ) return;
  used[id] = true;
  for(int i=0;i<(int)G[id].size();i++)
    dfs(G[id][i]);
  vs.push_back( id );
}
void rdfs(int id,int k){
  if( used[id] ) return;
  used[id] = true;
  cmp[id]=k;
  for(int i=0;i<(int)rG[id].size();i++)
    rdfs(rG[id][i],k);  
}

int scc(int size){
  memset(used,0,sizeof(used));
  for(int i=0;i<size;i++)
    dfs(i);
  memset(used,0,sizeof(used));
  int k=0;
  for(int i=(int)vs.size()-1;i>-1;i--){
    if( used[vs[i]] ) continue;
    rdfs(vs[i],k++);
  }  
  return k;
}
//ここまで------------------------------------------------------

//a->~aを返す,ただし頂点はa,~a,b,~b,c,~c...,の順に0,1,2,3,4,5...と割り振る
int rev(int a){
  if( a & 1 ) return --a;
  else return ++a;
}

//2-SAT本体
//頂点はa,~a,b,~b,c,~c...,の順に0,1,2,3,4,5...と割り振る
//( v[0].first or v[0].second ) and ( v[1].first or v[i].second ) and ...
bool two_SAT(const vector<P>& v,int size){
  init( size );
  for(int i=0;i<(int)v.size();i++){
    //aVb
    add_edge( rev(v[i].first), v[i].second );//~a->b
    add_edge( rev(v[i].second), v[i].first );//~b->a
  }
  scc(size);
  //xと~xが異なる強連結成分に含まれるか判定
  for(int i=0;i<size;i+=2){
    if( cmp[i] == cmp[i+1] ) return false;
  }
  return true;
}