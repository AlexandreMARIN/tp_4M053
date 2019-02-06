#include "include/iterative_solvers.hpp"

using namespace std;

int main(){


  cout << "We solve the linear system Ax = b where:\nA :=\n";
  Matrix A(5);
  A.diag(2.);
  A.diag(-1.0, 1);
  A.diag(-1.0, -1);
  A.display();

  Vect b({1.0, 1.0, 1.0, 1.0, 1.0});
  cout << "\nb :=\n";
  b.display();

  Jacobi_.set_tol(1e-5);
  Jacobi_.set_A(&A);
  Jacobi_.set_b(&b);
  Jacobi_.solve();
  cout << "with Jacobi's method:\n\niterations : " << Jacobi_.get_niter() << "\n\nx = \n";
  Jacobi_.get_x().display();


  GaussSeidel_.set_tol(1e-5);
  GaussSeidel_.set_A(&A);
  GaussSeidel_.set_b(&b);
  GaussSeidel_.solve();
  cout << "with Gauss-Seidel's method:\n\niterations : " << GaussSeidel_.get_niter() << "\n\nx = \n";
  GaussSeidel_.get_x().display();


  Relax_.set_tol(1e-5);
  Relax_.set_A(&A);
  Relax_.set_b(&b);
  Relax_.set_omega(1.4);
  Relax_.solve();
  cout << "with Relaxation's method:\n\niterations : " << Relax_.get_niter() << "\n\nx = \n";
  Relax_.get_x().display();

  return 0;
}
