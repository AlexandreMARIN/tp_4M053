#include "include/Vect.hpp"

using namespace std;

int main(){

  Vect v;
  Vect v2(5);
  Vect v3(v2);

  cout << v;

  v2.display();

  v3(0) = 5.0;

  cout << 2.0*v3;


  return 0;
}
