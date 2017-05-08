vector<int> G[MAX_V];

int pa[MAX_V],de[MAX_V],sz[MAX_V];
int cd[MAX_V],mv[MAX_V],hl[MAX_V];

void HLD(int root=0){
  queue<int> que;

  for(int i=1;i<=V;i++){
    sz[i]=1;
    pa[i]=cd[i]=mv[i]=i;
    hl[i]=G[i].size();
    if(i!=root&&hl[i]==1)que.push(i);
  }
  // calc sz[] cd[] pa[]
  while(!que.empty()){
    int pos=que.front();que.pop();
    if(pos==root)continue;
    for(int i=0;i<(int)G[pos].size();i++){
      int to=G[pos][i];
      if(pa[to]==pos)continue;
      pa[pos]=to;
      sz[to]+=sz[pos];           
      if(cd[to]==to||sz[cd[to]]<sz[pos])cd[to]=pos;
      hl[to]--;
      if(hl[to]==1)que.push(to);
    }
  }

  // calc hl[] de[] mv[] 
  que.push(root);
  de[root]=hl[root]=0;
  while(!que.empty()){
    int pos=que.front();que.pop();
    if(cd[pos]==pos)continue;
    mv[cd[pos]]=mv[pos];
    hl[cd[pos]]=hl[pos]+1;
    int sum=hl[pos]+1+sz[cd[pos]];
    for(int i=0;i<(int)G[pos].size();i++){
      int to=G[pos][i];
      if(to==pa[pos])continue;
      que.push(to);
      de[to]=de[pos]+1;
      if(to==cd[pos])continue;
      hl[to]=sum;
      sum+=sz[to];
    }
  }
}
  
int lca(int a,int b){
  while(mv[a]!=mv[b]){
    if(de[mv[a]]<de[mv[b]])swap(a,b);
    a=pa[mv[a]];
  }
  return (de[a]<de[b]?a:b);
}

void func(edge e){
  int x=e.from,y=e.to;
  while(mv[x]!=mv[y]){
    if(de[mv[x]]<de[mv[y]])swap(x,y);
    T.set(hl[mv[x]],hl[x]+1,e.cost);
    x=pa[mv[x]];
  }
  if(hl[x]>hl[y])swap(x,y);
  if(hl[x]+1<hl[y]+1){
    T.set(hl[x]+1,hl[y]+1,e.cost);
  }
}
