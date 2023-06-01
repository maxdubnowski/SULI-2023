#include <iostream>
#include <string>
#include <sstream>
using namespace std;

struct movies{
  string title;
  int cost;
  int year;
} mine, yours;

void printMovie(movies movie){
  cout << movie.title << " ("<< movie.year << ") - $" <<movie.cost << endl;;
}

int main (){
  string mystr;
  mine.title = "The Joker";
  mine.year = 2018;
  mine.cost = 12;

  cout << "Enter your Favorite Title: " ;
  getline (cin, yours.title);
  cout << "Enter year: ";
  getline(cin, mystr);
  stringstream(mystr) >> yours.year;
  
  cout << "Enter cost: ";
  cin >> mystr;
  stringstream(mystr) >> yours.cost;

  cout << "My favorite movie is: " <<endl;
  printMovie(mine);

  cout << "Yours is: " <<endl;
  printMovie(yours);

  movies *pmovie;
  pmovie = &mine;

  cout << "Pointer Movie title: " << pmovie->title << endl;
  cout << "Pointer Movie cost: " << (*pmovie).title << endl; //equivalent to above method

  return 0;
}

