#include <fstream>
#include <cmath>

#include "../include/iterative_solvers.hpp"

using namespace std;

Matrix create_A_N(int N){
  Matrix A_N(N);

  A_N.diag(2.0);
  A_N.diag(-1.0, 1);
  A_N.diag(-1.0, -1);

  return A_N;
}

int main(){

  vector<int> sizes({2, 4, 8, 16, 32, 64, 128});

  Jacobi_.set_tol(0.1);
  Jacobi_.set_n_max(100000);
  ofstream file("results.data");
  Matrix A(0);
  Vect b(0);



  Jacobi_.set_A(&A);
  Jacobi_.set_b(&b);


  for(unsigned int i=0;i<sizes.size();i++){
    A = create_A_N(sizes[i]);
    b.resize(sizes[i]);
    b.fill(1.0);

    file << sizes[i] << " ";
    Jacobi_.solve();
    file << Jacobi_.get_niter() << " ";

    //"analytical assessment"
    file << (log(0.1)/log(cos(M_PI/(sizes[i]+1)))) << " ";

    //assessment no 3
    file << -2.0*log(0.1)*(sizes[i]+1)*(sizes[i]+1)/(M_PI*M_PI) << "\n";
  }





  return 0;
}
