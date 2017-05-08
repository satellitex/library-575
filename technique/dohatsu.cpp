#include<bits/stdc++.h>
using namespace std;
typedef long long ll;

ll largep   = 105757705751LL; // 12 digit
ll prime575 =  57557557999LL; // 11 digit
ll prime579 =   5755755799LL; // 10 digit



ll smallprime = 1234567891; // 10 digit


//O(1), ただし、m <= 2^42 (4.39*10^12)


ll modmul(ll a, ll b, ll m){
  return (((a * (b >> 20) % m) << 20) + a * (b & ((1 << 20) - 1))) % m;
}

bool isP(ll x){
  for(ll i=2;i*i<=x;i++){
    if(x%i==0)return false;
  }
  return true;
}

int main(){
  cout<< isP(smallprime) <<endl;
  cout<< smallprime <<endl;
  cout<< smallprime*smallprime <<endl;
  return 0;
}
