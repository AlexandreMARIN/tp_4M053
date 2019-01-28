#include "include/Matrix.hpp"

using namespace std;


int main(){

  Matrix A, B(3);
  Matrix C(A);
  Vect v(3);

  v(0) = 1.;
  v(1) = 2.;
  v(2) = -1.;

  cout << "size of A : " << A.size() << "\n";

  A.display();

  A(0, 0) = 1.;
  A(1, 1) = 2.;
  A(2, 2) = 3.;
  B(0, 0) = 1.5;
  B(1, 1) = 2.;
  B(2, 2) = .5;

  cout << "A = \n" << A << "B = " << B;

  cout << "A+B = \n" << A+B;
  C = A - B;
  cout << "A-B = \n";
  C.display();
  cout << "A*B = \n" << A*B;
  cout << "2.*A = \n" << 2.*A;

  cout << "v : " << v;
  cout << "Av = \n" << A*v;

  return 0;
}
