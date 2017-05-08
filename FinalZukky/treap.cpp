#include<bits/stdc++.h>
using namespace std;

#define INF (1000000009)
#define MAX (1<<19)

typedef unsigned int uint;
extern inline uint xor128(void){
  static uint x=123456789,y=362436069,z=521288629,w=88675123; 
  uint t=x^(x<<11);
  x=y; y=z; z=w;
  return w=(w^(w>>19))^(t^(t>>8)); 
}

//�m�[�h
struct node {
  node *l,*r;
  int val;//�l
  int pri;//�D��x
  int cnt;//�����؂̃T�C�Y
  int sum;//�����؂̒l�̘a
  node(int v) : val(v), cnt(1), sum(v) {
    l = r = NULL;
    pri = xor128();
  }
};

typedef pair<node*,node*> pnn;
//���s�񕪒T����
struct Treap{
  node* root;
  Treap(){ root = NULL; }

  int count( node * t){ return !t ? 0:t->cnt; }
  int sum( node *t ) { return !t ? 0:t->sum; }

  node *update(node *t ){
    t->cnt = count(t->l) + count(t->r) + 1;
    t->sum = sum(t->l) + sum(t->r) + 1;
    return t;
  }

  node *merge(node *l, node *r ){
    if( !l || !r ) return !l ? r : l;

    if( l->pri > r ->pri ){
      l->r = merge(l->r,r);
      return update(l);
    } else {
      r->l = merge(l,r->l);
      return update(r);
    }
  }

  pnn split(node* t,int k){
    if( !t ) return pnn(NULL,NULL);
    if( k <= count(t->l) ){
      pnn s = split(t->l,k);
      t->l = s.second;
      return make_pair(s.first,update(t));
    } else {
      pnn s = split(t->r,k - count(t->l) - 1);
      t->r = s.first;
      return make_pair(update(t),s.second);
    }
  }

  int search(int v){//v�ȉ����������邩
    node* t = root;
    int sum = 0;
    int id = 0;
    while( t ){
      id = sum + count(t->l);
      if( t->val <= v ){
	sum += count(t->l) + 1;
	t = t->r;
      } else {
	t = t->l;
      }
    }
    return sum;
  }

  bool isexist(int v){//v�����邩�ۂ�
    node *t = root;
    while( t ){
      if( t->val < v ){
	t = t->r;
      } else if( t->val > v ){
	t = t->l;
      } else
	return true;
    }
    return false;
  }

  node* insert(node *t, int k, int v ){
    pnn s = split(t,k);
    s.first = merge(s.first,new node(v));
    return merge(s.first,s.second);
  }

  node* insert(int v){
    if( !root ) return root = new node(v);
    return root = insert( root, search(v), v );
  }

  void print(node* v){
    if( !v ) return;
    print(v->l);
    cout << v->val << " ";
    print(v->r);
  }
  void print(){
    cout << "print :";
    print(root);
    cout << endl;
  }
  
};
