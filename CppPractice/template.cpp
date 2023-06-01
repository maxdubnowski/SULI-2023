#include <iostream>
using namespace std;

template <class T>
T GetMax (T a, T b){
  return ((a>b)? a:b);
}

int main(){
  int i=5, j=6, k;
  float m= 3.14, n =2.71, p;

  k = GetMax<int>(i,j);
  p = GetMax<float>(m,n);

  cout << k << endl;
  cout << p << endl;


  return 0;
}
