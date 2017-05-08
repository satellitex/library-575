#include<bits/stdc++.h>
using namespace std;
#define INF (1<<30)

//�x���]���Z�O�����g�c���[�i�����ԂɈ�l�ɑ����A��l��0�ŏ������@�{�@��Ԃ̍ő�l or �ŏ��l or ���a �����߂�@���삪�\)
//��l��0�ɏ������Ɋւ��Ă͎኱O���d���Ȃ�C������̂Ŏg��Ȃ����͏����Ă��������B
//�w�肷���Ԃ͏�� [a,b) �̔��J��Ԃł��邱�Ƃɒ��ӂ��ĉ������B
struct segtree{
  //�ő�l�����߂�p
  vector<int> datamax;

  //�ŏ��l�����߂�p
  vector<int> datamin;

  //���a�����߂�p
  vector<int> datasum;

  //�x���p
  vector<int> delay;

  int n;

  //������
  void init(int _n){
    n = 1;
    while( n < _n ) n*=2;
    datamax.resize( 2 * n );
    datamin.resize( 2 * n );
    datasum.resize( 2 * n );
    delay.resize( 2 * n );
  }

  //�N���A(���x���g���ꍇ�j
  void clear(){
    datamax.clear();
    datamin.clear();
    datasum.clear();
    delay.clear();
  }
 
  //	 �x���]���{��
  /*	l,r �Ɋւ��Ă͑��a�����߂鎞�ɂ̂ݎg���Ă���B */
  /*	0�ŏ��������铮�삪����Ȃ��ꍇ�͎O�����Z�q�̕�����
  //    if( delay[2*k+i]<0 && delay[k] > 0 ) delaycalc(2*k+i, i==0?l:(l+r)/2, i==0?(l+r)/2:r) �̕���������Ȃ�
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

  //�ȈՔŁi���a�ɑ΂��鏈������and0����������ver)
  void delaycalc(int k){
    datamax[k] = datamax[k]+delay[k];
    datamin[k] = datamin[k]+delay[k];
    if( k+1 < n )
      for(int i=1;i<=2;i++)
	delay[2*k+i] = delay[2*k+i]+delay[k];
    delay[k] = 0;
  }
  

  //��ԁ@[a,b)�@�Ɉ�l�Ɂ@x�@�𑫂� 
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

  //�����Ԃ̍ő�l�����߂�
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

  //�����Ԃ̍ŏ��l�����߂�
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

  //�����Ԃ̑��a�����߂�
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

�@//�e�ȈՃN�G���Ăяo��
  void reset(int a,int b){ add( a, b, -1,  0, 0, n ); }//0�ŏ�����
  void add( int a,int b,int x){ add( a, b, x, 0,0,n); }//������
  int querymax(int a,int b){ return querymax(a,b,0,0,n); }//�ő�l�擾
  int querymin(int a,int b){ return querymin(a,b,0,0,n); }//�ŏ��l�擾
  int querysum(int a,int b){ return querysum(a,b,0,0,n); }//���a�擾
};
