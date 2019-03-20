#include "./include/conj_grad.hpp"

#include <iostream>

using namespace std;

int main(){

  int n, n_max;

  cout << "Give an integer : ";
  cin >> n;

  if(cin.bad()){
    cout << "\nn = 2\n";
    n = 2;
  }

  cout << "Number of iterations : ";
  cin >> n_max;

  if(cin.bad()){
    cout << "\nn_max = \n" << n;
    n_max = n;
  }

  Matrix A(n);
  A.diag(2.0);
  A.diag(-1.0, 1);
  A.diag(-1.0, -1);

  Vect b(n);
  b.fill(1.0);


  ConjGrad_.set_A(&A);
  ConjGrad_.set_b(&b);
  ConjGrad_.set_tol(1e-2);
  ConjGrad_.set_n_max(n_max);
  ConjGrad_.solve();

  cout << "The linear system Ax=b with\nA:=\n";
  A.display();
  cout << "and b:=\n";
  b.display();
  cout << "has for solution x = \n";
  ConjGrad_.get_x().display();

  cout << "solution found with " << ConjGrad_.get_niter() << " iterations\n";

  return 0;
}
