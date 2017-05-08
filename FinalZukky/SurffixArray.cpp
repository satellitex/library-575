#include <string>
#include <vector>
#include <algorithm>
using namespace std;

const int MAX_N = 100000;

//ここから接尾辞配列(suffix_array)関連（classにするとcompareでエラーが出たため、グローバル変数に）
int n,k;
int rank[MAX_N+1];
int tmp[MAX_N+1];

bool compare_sa(int i, int j){
  if(rank[i] != rank[j]) return rank[i] < rank[j];
  else{
    int ri = i + k <= n ? rank[i+k] : -1;
    int rj = j + k <= n ? rank[j+k] : -1;
    return ri < rj;
  }
}

/*
  文字列strの接尾辞配列(sa)の構築
  O(n * log^2(n))
*/
void construct_sa(string str, int *sa){
  n = str.length();
  
  for(int i=0;i<=n;i++){
    sa[i] = i;
    rank[i] = i < n ? str[i] : -1;
  }
  
  for(k=1;k<=n;k*=2){
    sort(sa, sa+n+1, compare_sa);
    tmp[sa[0]] = 0;
    for(int i=1;i<=n;i++) tmp[sa[i]] = tmp[sa[i-1]] + compare_sa(sa[i-1],sa[i]);
    for(int i=0;i<=n;i++) rank[i] = tmp[i];
  }
}

/*
  高さ配列(LCP)の構築
  sa[i]とsa[i+1]の文字列の先頭何文字が共通しているか（sa[0]は空）
  rankを再利用、O(n)
*/
void construct_lcp(string str, int *sa, int *lcp){
  for(int i=0;i<=n;i++) rank[sa[i]] = i;
  
  int h = 0;
  lcp[0] = 0;
  for(int i=0;i<n;i++){
    int j = sa[rank[i] - 1];
    if(h > 0) h--;
    for(; j + h < n && i + h < n; h++)
      if(str[j+h] != str[i+h]) break;
    lcp[rank[i] - 1] = h;
  }
}
