#include "include/iterative_solvers.hpp"

using namespace std;

int main(){

  Matrix A(5);
  A.diag(2.);
  A.diag(-1.0, 1);
  A.diag(-1.0, -1);
  A.display();
  Vect b({1.0, 1.0, 1.0, 1.0, 1.0});
  Jacobi_.set_A(&A);
  Jacobi_.set_b(&b);
  Jacobi_.solve();

  Jacobi_.get_x().display();

  return 0;
}
