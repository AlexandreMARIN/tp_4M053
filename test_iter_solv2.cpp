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
  double tol = 0.1;

  Jacobi_.set_tol(tol);
  Jacobi_.set_n_max(100000);

  GaussSeidel_.set_tol(tol);
  GaussSeidel_.set_n_max(100000);

  Relax_.set_tol(tol);
  Relax_.set_n_max(100000);



  ofstream file("Jacobi.data");
  Matrix A(0);
  Vect b(0);



  Jacobi_.set_A(&A);
  Jacobi_.set_b(&b);

  GaussSeidel_.set_A(&A);
  GaussSeidel_.set_b(&b);

  Relax_.set_A(&A);
  Relax_.set_b(&b);

  //Jacobi
  for(unsigned int i=0;i<sizes.size();i++){
    A = create_A_N(sizes[i]);
    b.resize(sizes[i]);
    b.fill(1.0);

    file << sizes[i] << " ";
    Jacobi_.solve();
    file << Jacobi_.get_niter() << " ";

    //"analytical assessment"
    file << (log(tol)/log(cos(M_PI/(sizes[i]+1)))) << " ";

    //assessment no 3
    file << -2.0*log(tol)*(sizes[i]+1)*(sizes[i]+1)/(M_PI*M_PI) << "\n";
  }

  file.close();

  file.open("GaussSeidel.data");
  //Gauss-Seidel
  for(unsigned int i=0;i<sizes.size();i++){
    A = create_A_N(sizes[i]);
    b.resize(sizes[i]);
    b.fill(1.0);

    file << sizes[i] << " ";
    GaussSeidel_.solve();
    cout << "A_*x_ = " << A*GaussSeidel_.get_x() << "\n\n";
    file << GaussSeidel_.get_niter() << " ";

    //"analytical assessment"
    file << (log(tol)/(2.0*log(cos(M_PI/(sizes[i]+1))))) << " ";

    //assessment no 3
    file << -log(tol)*(sizes[i]+1)*(sizes[i]+1)/(M_PI*M_PI) << "\n";
  }

  file.close();


  file.open("Relax.data");
  //Relax
  double omg_star, rho;
  for(unsigned int i=0;i<sizes.size();i++){
    A = create_A_N(sizes[i]);
    b.resize(sizes[i]);
    b.fill(1.0);
    rho = cos(M_PI/(sizes[i]+1));
    omg_star = 2.0/(1.0+sqrt( 1.0 - rho*rho ));
    cout << "omg_star(" << sizes[i] << ") = " << omg_star << "\n";
    Relax_.set_omega(omg_star);

    file << sizes[i] << " ";
    Relax_.solve();
    cout << "A_*x_ = " << A*Relax_.get_x() << "\n\n";
    file << Relax_.get_niter() << " ";

    //"analytical assessment"
    file << log(tol)/log(omg_star - 1.0) << " ";

    //assessment no 3
    file << -log(tol)*(sizes[i]+1)/(M_PI*M_PI) << "\n";
  }

  file.close();


  return 0;
}
