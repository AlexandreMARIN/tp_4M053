#include <fstream>
#include <cmath>


#include "../include/iterative_solvers.hpp"

using namespace std;
using namespace std::chrono;

Matrix create_A_N(int N){
  Matrix A_N(N);

  A_N.diag(2.0);
  A_N.diag(-1.0, 1);
  A_N.diag(-1.0, -1);

  return A_N;
}

int main(){

  ofstream file;
  constexpr int measures = 10;

  vector<int> sizes(measures);
  sizes[0] = 10;
  for(int i=1;i<measures;i++){
    sizes[i] = sizes[i-1]+10;
  }

  double tol = 1e-1;
  int n_max = 20000;
  Matrix A(0);
  Vect b(0);

  double omg_star, rho;

  Jacobi_.set_tol(tol);
  Jacobi_.set_n_max(n_max);
  Jacobi_.set_A(&A);
  Jacobi_.set_b(&b);

  GaussSeidel_.set_tol(tol);
  GaussSeidel_.set_n_max(n_max);
  GaussSeidel_.set_A(&A);
  GaussSeidel_.set_b(&b);

  Relax_.set_tol(tol);
  Relax_.set_n_max(n_max);
  Relax_.set_A(&A);
  Relax_.set_b(&b);



  //resvec
  A = create_A_N(200);
  b.resize(200);
  b.fill(1.0);
  vector<double> resvec;
  file.open("resvec.json");
  file << "{\n\t\"Jacobi\" : {\n\t\t\"resvec\" : [";
  Jacobi_.solve();
  resvec = Jacobi_.get_resvec();
  file << resvec[0];
  for(unsigned int i=1;i<resvec.size();i++){
    file << ", " << resvec[i];
  }
  file << "]\n\t},\n";

  file << "\n\t\"GaussSeidel\" : {\n\t\t\"resvec\" : [";
  GaussSeidel_.solve();
  resvec = GaussSeidel_.get_resvec();
  file << resvec[0];
  for(unsigned int i=1;i<resvec.size();i++){
    file << ", " << resvec[i];
  }
  file << "]\n\t},\n";

  file << "\n\t\"Relax\" : {\n\t\t\"resvec\" : [";
  rho = cos(M_PI/(201));
  omg_star = 2.0/(1.0+sqrt( 1.0 - rho*rho ));
  Relax_.set_omega(omg_star);
  Relax_.solve();
  resvec = Relax_.get_resvec();
  file << resvec[0];
  for(unsigned int i=1;i<resvec.size();i++){
    file << ", " << resvec[i];
  }
  file << "]\n\t}\n";
  file << "}";
  file.close();

  //time

  //Jacobi
  file.open("J_time");
  for(int s : sizes){
    A = create_A_N(s);
    b.resize(s);
    b.fill(1.0);

    file << s << " ";
    Jacobi_.solve();
    file << Jacobi_.get_duration() << "\n";

  }
  file.close();

  //GaussSeidel
  file.open("GS_time");
  for(int s : sizes){
    A = create_A_N(s);
    b.resize(s);
    b.fill(1.0);

    file << s << " ";
    GaussSeidel_.solve();
    file << GaussSeidel_.get_duration() << "\n";

  }
  file.close();

  //Relax
  file.open("R_time");

  for(int s : sizes){
    A = create_A_N(s);
    b.resize(s);
    b.fill(1.0);

    file << s << " ";
    rho = cos(M_PI/(s+1));
    omg_star = 2.0/(1.0+sqrt( 1.0 - rho*rho ));
    Relax_.set_omega(omg_star);
    Relax_.solve();
    file << Relax_.get_duration() << "\n";

  }
  file.close();


  return 0;
}
