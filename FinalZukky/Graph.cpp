#include <iostream>
#include <algorithm>
#include <vector>
#include <map>
#include <queue>
#include <functional>
using namespace std;

const int INF = 1<<28;

struct Edge {
  int to, cost, cap, rev;
  Edge(int to, int cost, int cap, int rev)
    : to(to), cost(cost), cap(cap), rev(rev) {}
  Edge() {}
};

typedef vector<vector<Edge> > Graph;

void addEdge(int from, int to, Graph &g) {
  g[from].push_back(Edge(to,0,0,-1));
}

void addEdge(int from, int to, int cost, Graph &g) {
  g[from].push_back(Edge(to,cost,0,-1));
}

void addEdgeF(int from, int to, int cost, int cap, Graph &g) {
  g[from].push_back(Edge(to,cost,cap,g[to].size()));
  g[to].push_back(Edge(from,-cost,0,(int)g[from].size()-1));
}

int dfsF(int v, int t, int f, vector<int> &used, Graph &g) {
  if(v == t) return f;
  used[v] = true;
  for(int i = 0; i < g[v].size(); ++i) {
    Edge &e = g[v][i];
    if(!used[e.to] && e.cap > 0) {
      int d = dfsF(e.to, t, min(f, e.cap), used, g);
      if(d > 0) {
        e.cap -= d;
        g[e.to][e.rev].cap += d;
        return d;
      }
    }
  }
  return 0;
}

int maxFlow(int s, int t, Graph g) {
  int flow = 0;
  while(1) {
    vector<int> used(g.size());
    int f = dfsF(s, t, INF, used, g);
    if(f == 0) break;
    flow += f;
  }
  return flow;
}

//こっから最大流の辺の削除と追加各O(E)-------------------------------------
int once_flow(int s,int t){
  memset(used,0,sizeof(used));
  return dfs(s,t,INF);
}
int N,V;
int S,T;
//辺の削除
int remove_flow(edge &u,edge &v){//u:正辺,v:逆辺
  bool f = false;
  if( u.cap == 0 ) f = true;
  u.cap = 0;
  v.cap = 0;
  if( f ){
    int fw = once_flow(v.to,u.to);
    if( fw == 0 ) {
      once_flow(T,u.to);    
      once_flow(v.to,S);
      if( once_flow(S,T) == 1 ) return 0;
      return -1;
    }
    return 0;
  }
  return 0;
}
//辺を0->1にする(いくつでもいいはず)
int add_flow(edge &u,edge &v){
  u.cap = 1;
  v.cap = 0;
  int ret =  once_flow(S,T);
  return ret;
}
//最大流押し戻しここまで----------------------------------



// O(EV^2)
namespace Dinic {
  Graph G;
  vector<int> level, iter;

  void bfs(int s) {
    level = vector<int>(G.size(), -1);
    queue<int> que;
    level[s] = 0;
    que.push(s);
    while(!que.empty()) {
      int v = que.front(); que.pop();
      for(int i = 0; i < G[v].size(); ++i) {
        const Edge &e = G[v][i];
        if(e.cap > 0 && level[e.to] < 0) {
          level[e.to] = level[v] + 1;
          que.push(e.to);
        }
      }
    }
  }

  int dfs(int v, int t, int f) {
    if(v == t) return f;
    for(int &i = iter[v]; i < G[v].size(); ++i) {
      Edge &e = G[v][i];
      if(e.cap > 0 && level[v] < level[e.to]) {
        int d = dfs(e.to, t, min(f, e.cap));
        if(d > 0) {
          e.cap -= d;
          G[e.to][e.rev].cap += d;
          return d;
        }
      }
    }
    return 0;
  }

  int maxFlow(int s, int t) {
    int flow = 0;
    while(1) {
      bfs(s);
      if(level[t] < 0) return flow;
      iter = vector<int>(G.size(), 0);
      for(int f; (f = dfs(s, t, INF)) > 0; ) {
        flow += f;
      }
    }
    return flow;
  }
}

typedef pair<int,int> Pair;

// O(FV^2)
int minCostFlow(int s, int t, int f, Graph g) {
  int V = g.size();
  vector<int> h(V), dist(V), prevv(V), preve(V);
  int res = 0;
  fill(h.begin(), h.end(), 0);
  while(f > 0) {
    priority_queue<Pair, vector<Pair>, greater<Pair> > que;
    fill(dist.begin(), dist.end(), INF);
    dist[s] = 0;
    que.push(Pair(0, s));
    while(!que.empty()) {
      Pair p = que.top(); que.pop();
      int v = p.second;
      if(dist[v] < p.first) continue;
      for(int i = 0; i < g[v].size(); ++i) {
        Edge &e = g[v][i];
        int tmp = dist[v] + e.cost + h[v] - h[e.to];
        if(e.cap > 0 && dist[e.to] > tmp) {
          dist[e.to] = tmp;
          prevv[e.to] = v;
          preve[e.to] = i;
          que.push(Pair(dist[e.to], e.to));
        }
      }
    }
    if(dist[t] == INF) {
      return -1;
    }
    for(int v = 0; v < V; ++v) h[v] += dist[v];
    
    int d = f;
    for(int v = t; v != s; v = prevv[v]) {
      d = min(d, g[prevv[v]][preve[v]].cap);
    }
    f -= d;
    res += d * h[t];
    for(int v = t; v != s; v = prevv[v]) {
      Edge &e = g[prevv[v]][preve[v]];
      e.cap -= d;
      g[v][e.rev].cap += d;
    }
  }
  return res;
}

// from topologicalSort
bool dfsTS(const Graph &g, int v, int &times, vector<int> &vis, 
           vector<int> &num, vector<int> &res) {
  vis[v] = true;
  for(int i = 0; i < g[v].size(); ++i) {
    int nv = g[v][i].to;
    if(vis[nv]) {
      if(!num[nv]) return false;
      continue;
    }
    if(!dfsTS(g, nv, times, vis, num, res)) return false;
  }
  num[v] = ++times;
  res.push_back(v);
  return true;
}

// トポロジカルソートされた頂点のリストを返す。
// 循環がある場合、空のvectorを返す。
// O(V+E)
// verified by 
vector<int> topologicalSort(const Graph &g) {
  int n = g.size();
  vector<int> vis(n), num(n), res;
  int times = 0;
  for(int i = 0; i < n; ++i) {
    if(vis[i]) continue;
    if(!dfsTS(g, i, times, vis, num, res)) return vector<int>();
  }
  reverse(res.begin(), res.end());
  return res;
}

namespace LCA {
  // input
  int V;
  Graph G;
  int root;

  const int MAX_LOG_V = 32;
  vector<vector<int> > parent;
  vector<int> depth;

  void dfs(int v, int p, int d) {
    parent[0][v] = p;
    depth[v] = d;
    for(int i = 0; i < G[v].size(); ++i) {
      if(G[v][i].to != p) dfs(G[v][i].to, v, d + 1);
    }
  }

  void init() {
    parent = vector<vector<int> >(MAX_LOG_V, vector<int>(V));
    depth  = vector<int>(V);

    dfs(root, -1, 0);
    for(int k = 0; k + 1 < MAX_LOG_V; ++k) {
      for(int v = 0; v < V; ++v) {
        if(parent[k][v] < 0) parent[k+1][v] = -1;
        else parent[k+1][v] = parent[k][parent[k][v]];
      }
    }
  }

  int lca(int u, int v) {
    if(depth[u] > depth[v]) swap(u, v);
    for(int k = 0; k < MAX_LOG_V; ++k) {
      if((depth[v] - depth[u]) >> k & 1) {
        v = parent[k][v];
      }
    }
    if(u == v) return u;
    for(int k = MAX_LOG_V - 1; k >= 0; --k) {
      if(parent[k][u] != parent[k][v]) {
        u = parent[k][u];
        v = parent[k][v];
      }
    }
    return parent[0][u];
  }
}

// from getBridge
int dfsB(const Graph &g, int v, int prev, int &times, vector<int> &num,
         vector<pair<int,int> > &bridge) {
  int res = num[v] = times++;
  for(int i = 0; i < g[v].size(); ++i) {
    int nv = g[v][i].to;
    int t;
    if(num[nv] == -1) {
      t = dfsB(g, nv, v, times, num, bridge);
      if(t > num[v]) {
        bridge.push_back(make_pair(v, nv));
      }
    } else {
      t = num[nv];
    }
    if(nv != prev) {
      res = min(res, t);
    }
  }
  return res;
}

// 橋を求める。
// O(V+E)
// 多重辺があった場合はそれを１本の辺として見るので注意。
// verified by UVa796 2013/09/28
//
// また、元のグラフから橋を取り除いたグラフにある
// 各連結成分は二重辺連結成分である。
// verified by LiveArchive 6044 2014/10/14
vector<pair<int,int> > getBridge(const Graph &g) {
  vector<pair<int,int> > res;
  int n = g.size();
  int times = 0;
  vector<int> num(n, -1);
  for(int i = 0; i < n; ++i) {
    if(num[i] == -1) {
      dfsB(g, i, -1, times, num, res);
    }
  }
  return res;
}

// from getArticulationPoint
int dfsA(const Graph &g, int v, int prev, int &times, vector<int> &num,
         vector<int> &ap) {
  int res = num[v] = times++;
  bool flag = false;
  int c = 0;
  for(int i = 0; i < g[v].size(); ++i) {
    int nv = g[v][i].to;
    int t;
    if(num[nv] == -1) {
      t = dfsA(g, nv, v, times, num, ap);
      ++c;
      if(~prev && t >= num[v]) {
        flag = true;
      }
    } else {
      t = num[nv];
    }
    if(nv != prev) {
      res = min(res, t);
    }
  }
  if(prev == -1 && c > 1) flag = true;
  if(flag) ap.push_back(v);
  return res;
}

// 関節点を求める。
// O(V+E)
// 多重辺があった場合はそれを１本の辺として見るので注意。
// verified by 
vector<int> getArticulationPoint(const Graph &g) {
  vector<int> res;
  int n = g.size();
  int times = 0;
  vector<int> num(n, -1);
  for(int i = 0; i < n; ++i) {
    if(num[i] == -1) {
      dfsA(g, i, -1, times, num, res);
    }
  }
  return res;
}

// from scc
void dfsC(int v, const Graph &g, vector<int> &used, vector<int> &vs) {
  used[v] = true;
  for(int i = 0; i < g[v].size(); ++i) {
    if(!used[g[v][i].to]) dfsC(g[v][i].to, g, used, vs);
  }
  vs.push_back(v);
}

// from scc
void rdfsC(int v, int k, const Graph &rg, vector<int>&used, vector<int> &cmp) {
  used[v] = true;
  cmp[v] = k;
  for(int i = 0; i < rg[v].size(); ++i) {
    if(!used[rg[v][i].to]) rdfsC(rg[v][i].to, k, rg, used, cmp);
  }
}

// 強連結成分(SCC:Strongly Connected Component)を求める。
// O(V+E)
// from Ari Book
vector<int> scc(const Graph &g) {
  int V = g.size();
  Graph rg(V);
  for(int i = 0; i < V; ++i) {
    for(int j = 0; j < g[i].size(); ++j) {
      addEdge(g[i][j].to, i, rg);
    }
  }
  vector<int> used(V, 0);
  vector<int> vs;
  for(int v = 0; v < V; ++v) {
    if(!used[v]) dfsC(v, g, used, vs);
  }
  used = vector<int>(V, 0);
  vector<int> cmp(V);
  int k = 0;
  for(int i = (int)vs.size() - 1; i >= 0; --i) {
    if(!used[vs[i]]) rdfsC(vs[i], k++, rg, used, cmp);
  }
  return cmp;
}


//2-SAT (sccを利用する)--------------------------------------------------------------
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
//ここまで-----------------------------------------------------------------------------




// uva 10330
int main() {
  int N, M, B, D;
  Graph g;
  while(cin >> N) {
    g = Graph(N*2+10);
    for(int i = 0; i < 500; ++i) g[i].clear();
    for(int i = 1; i <= N; ++i) {
      int cap;
      cin >> cap;
      addEdgeF(i*2, i*2+1, 0, cap, g);
    }
    cin >> M;
    for(int i = 0; i < M; ++i) {
      int from, to, cap;
      cin >> from >> to >> cap;
      addEdgeF(from*2+1, to*2, 0, cap, g);
    }
    cin >> B >> D;
    for(int i = 0; i < B; ++i) {
      int to;
      cin >> to;
      addEdgeF(0*2+1, to*2, 0, INF, g);
    }
    for(int i = 0; i < D; ++i) {
      int from;
      cin >> from;
      addEdgeF(from*2+1, (N+1)*2, 0, INF, g);
    }
    cout << maxFlow(0*2+1, (N+1)*2, g) << endl;
  }
  return 0;
}
