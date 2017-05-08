#include<bits/stdc++.h>
using namespace std;
#define INF (1<<30)

//遅延評価セグメントツリー（ある区間に一様に足す、一様に0で初期化　＋　区間の最大値 or 最小値 or 総和 を求める　動作が可能)
//一様に0に初期化に関しては若干Oが重くなる気がするので使わない時は消してください。
//指定する区間は常に [a,b) の半開区間であることに注意して下さい。
struct segtree{
  //最大値を求める用
  vector<int> datamax;

  //最小値を求める用
  vector<int> datamin;

  //総和を求める用
  vector<int> datasum;

  //遅延用
  vector<int> delay;

  int n;

  //初期化
  void init(int _n){
    n = 1;
    while( n < _n ) n*=2;
    datamax.resize( 2 * n );
    datamin.resize( 2 * n );
    datasum.resize( 2 * n );
    delay.resize( 2 * n );
  }

  //クリア(何度も使う場合）
  void clear(){
    datamax.clear();
    datamin.clear();
    datasum.clear();
    delay.clear();
  }
 
  //	 遅延評価本体
  /*	l,r に関しては総和を求める時にのみ使っている。 */
  /*	0で初期化する動作がいらない場合は三項演算子の部分と
  //    if( delay[2*k+i]<0 && delay[k] > 0 ) delaycalc(2*k+i, i==0?l:(l+r)/2, i==0?(l+r)/2:r) の部分がいらない
  //                                                  */
  void delaycalc(int k,int l,int r){
    datamax[k] = delay[k]<0?0:datamax[k]+delay[k];
    datamin[k] = delay[k]<0?0:datamin[k]+delay[k];
    datasum[k] = delay[k]<0?0:datasum[k] + delay[k]*(r-l);
    if( k+1 < n ){
      for(int i=1;i<=2;i++){
	if( delay[2*k+i]<0 && delay[k] > 0 ) delaycalc(2*k+i, i==1?l:(l+r)/2, i==1?(l+r)/2:r);
	delay[2*k+i] = delay[k]<0?-1:delay[2*k+i]+delay[k];
      }
    } 
    delay[k] = 0;
  }

  //簡易版（総和に対する処理無しand0初期化無しver)
  void delaycalc(int k){
    datamax[k] = datamax[k]+delay[k];
    datamin[k] = datamin[k]+delay[k];
    if( k+1 < n )
      for(int i=1;i<=2;i++)
	delay[2*k+i] = delay[2*k+i]+delay[k];
    delay[k] = 0;
  }
  

  //区間　[a,b)　に一様に　x　を足す 
  void add(int a,int b,int x,int k,int l,int r){
    delaycalc(k,l,r);
    if( r<=a || b<=l ) return;
    else if( a<=l && r<=b ){
      delay[k] = x;
      delaycalc(k,l,r);
    } else {
      add( a, b, x, 2*k + 1, l, (l+r)/2 );
      add( a, b, x, 2*k + 2, (l+r)/2, r );
      datamax[k] = max( datamax[2*k+1], datamax[2*k+2] );
      datamin[k] = min( datamin[2*k+1], datamin[2*k+2] );
      datasum[k] = datasum[2*k+1]+datasum[2*k+2];
    }
  }

  //ある区間の最大値を求める
  int querymax(int a,int b,int k,int l,int r){
    delaycalc(k,l,r);
    if( r<=a || b<=l ) return 0;
    if( a<=l && r<=b ){
      return datamax[k]+delay[k];
    } else {
      int vl = querymax( a, b, 2*k+1, l,(l+r)/2 );
      int vr = querymax( a, b, 2*k+2, (l+r)/2,r );
      return max(vl,vr);
    }
  }

  //ある区間の最小値を求める
  int querymin(int a,int b,int k,int l,int r){
    delaycalc(k,l,r);
    if( r<=a || b<=l ) return INF;
    if( a<=l && r<=b ){
      return datamin[k]+delay[k];
    } else {
      int vl = querymin( a, b, 2*k+1, l,(l+r)/2 );
      int vr = querymin( a, b, 2*k+2, (l+r)/2,r );
      return min(vl,vr);
    }
  }

  //ある区間の総和を求める
  int querysum(int a,int b,int k,int l,int r){
    delaycalc(k,l,r);
    if( r<=a || b<=l ) return 0;
    if( a<=l && r<=b ){
      return datasum[k]+delay[k]*(r-l);
    } else {
      int vl = querysum( a, b, 2*k+1, l,(l+r)/2 );
      int vr = querysum( a, b, 2*k+2, (l+r)/2,r );
      return vl+vr;
    }
  }

　//各簡易クエリ呼び出し
  void reset(int a,int b){ add( a, b, -1,  0, 0, n ); }//0で初期化
  void add( int a,int b,int x){ add( a, b, x, 0,0,n); }//加える
  int querymax(int a,int b){ return querymax(a,b,0,0,n); }//最大値取得
  int querymin(int a,int b){ return querymin(a,b,0,0,n); }//最小値取得
  int querysum(int a,int b){ return querysum(a,b,0,0,n); }//総和取得
};
