#include "include/Matrix.hpp"

using namespace std;

DenseMatrix create_A_N(int N){
  DenseMatrix A_N(N);

  A_N.diag(2.0);
  A_N.diag(-1.0, 1);
  A_N.diag(-1.0, -1);

  return A_N;
}


int main(){

  int n;

  cout << "Give an integer : ";
  cin >> n;
  if(cin.fail() || n==0){
    cout << "the integer 1 is chosen\n";
    n = 1;
  }

  DenseMatrix A(create_A_N(n));


  cout << "A :=\n";
  A.display();

  A.cholesky();

  A.display_Cholesky();

  Vect b(n);
  b.fill(1.0);

  cout << "the solution of the linear system Ax=b, with b = ";
  b.display();

  Vect x(A.solve_via_Cholesky(b));

  cout << "\n is x = \n";
  x.display();

  return 0;
}
