#include <bits/stdc++.h>
using namespace std;

// サイコロ関連
enum FACE { TOP, BOTTOM, FRONT, BACK, LEFT, RIGHT };

struct dice{
public:
  dice(){//ここ順不動,enumにあわせる
    id[0] = TOP;
    id[1] = FRONT;
    id[2] = LEFT;
    id[3] = RIGHT;
    id[4] = BACK;
    id[5] = BOTTOM;
  }
  
  int& operator[] (FACE f){return var[id[f]];}
  const int& operator[] (FACE f)const {return var[id[f]];}

  bool operator==(const dice &a)const {
    const dice &b = *this;
    for(int i=0;i<6;i++){
      if(a[(FACE)i] != b[(FACE)i]) return false;
    }
    return true;
  }

  void roll_x(){roll(TOP,BACK,BOTTOM,FRONT);} // X軸を基に回転 (FRONT -> BOTTOM)
  void roll_y(){roll(TOP,LEFT,BOTTOM,RIGHT);} // Y軸を基に回転 (TOP -> RIGHT)
  void roll_z(){roll(FRONT,RIGHT,BACK,LEFT);} // Z軸を基に回転 (FRONT -> LEFT)

  // 全通り（24通り）の回転
  vector<dice> all_rolls(){
    vector<dice> res;
    for(int k=0; k<6; (k&1 ? roll_y() : roll_x()), k++){
      for(int i=0; i<4; roll_z(),i++){
        res.push_back(*this);
      }
    }
    return res;
  }

  // 向き違いを考慮したequal
  bool equivalent_to(const dice &di){
    for(int k=0; k<6; (k&1 ? roll_y() : roll_x()), k++){
      for(int i=0; i<4; roll_z(),i++){
        if(*this == di)return true;
      }
    }
    return false;
  }

private:
  int var[6]; // 実際の値の記憶
  int id[6]; // 回転の記憶
  void roll(FACE a,FACE b,FACE c,FACE d){//a b c d →b c d a
    int tmp = id[a];
    id[a] = id[b];
    id[b] = id[c];
    id[c] = id[d];
    id[d] = tmp;
  }
};
