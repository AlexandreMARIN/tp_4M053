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

  A.decomp_LU();

  A.display_LU();

  Vect b(n);
  b.fill(1.0);

  cout << "the solution of the linear system Ax=b, with b = ";
  b.display();

  Vect x(A.solve_via_LU(b));

  cout << "\n is x = \n";
  x.display();


  cout << "-----------------------------";
  DenseMatrix B(n);
  B.diag(1.0);
  B(n-1, 0) = 2.0;
  B(0, n-1) = 3.0;
  B((n-1)/2 + 1, (n-1)/2) = 5.0;

  cout << "\nB :=\n";
  B.display();

  B.decomp_LU();

  B.display_LU();

  return 0;
}
