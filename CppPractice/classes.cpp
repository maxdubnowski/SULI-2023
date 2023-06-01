#include <iostream>
using namespace std;

class Rectangle{
  int width, height; //private

public:
  Rectangle(){ //default constructor
    width = 5;
    height =5;
  }
  
  Rectangle(int a, int b) { //overloading constructors 
    width = a;
    height = b;
  }

  int Area(void){ //defining a function inside the class
    return width*height;
  }

  int Perimeter (void);
};
  
int Rectangle::Perimeter(void){ //going into the scope of the Rectangle class
  return 2*width+2*height;
}



int main(){
  Rectangle recA(5,3);
  Rectangle recB;
  cout << "recA Area: " << recA.Area() << endl;
  cout << "recB Area: " << recB.Area() << endl;

  return 0;
}

