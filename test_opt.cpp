#include "./include/opt_step_grad.hpp"

#include <iostream>

using namespace std;

int main(){

  int n;

  cout << "Give an integer : ";
  cin >> n;

  if(cin.bad()){
    cout << "\nn = 2\n";
    n = 2;
  }

  DenseMatrix A(n);
  A.diag(2.0);
  A.diag(-1.0, 1);
  A.diag(-1.0, -1);

  Vect b(n);
  b.fill(1.0);


  OptStepGrad_.set_A(&A);
  OptStepGrad_.set_b(&b);
  OptStepGrad_.solve();

  cout << "The linear system Ax=b with\nA:=\n";
  A.display();
  cout << "and b:=\n";
  b.display();
  cout << "has for solution x = \n";
  OptStepGrad_.get_x().display();

  cout << "solution found with " << OptStepGrad_.get_niter() << " iterations\n";

  return 0;
}
