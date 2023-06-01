#include <stdio.h>
#include <iostream>
using namespace std;

int main(){
  int numbers[5];
  int *p;
  p = numbers; *p = 10;
  p = &numbers[1]; *p= 20;
  p++; *p=30;
  p = numbers +3; *p = 40;
  p = numbers; *(p+4) = 50;

  for (int i=0; i<5; i++){
    cout << numbers[i] << endl;
  }

  cout << "Go Bananas" << endl;
  return 0;
}
