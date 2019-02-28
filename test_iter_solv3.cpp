#include <fstream>
#include <cmath>


#include "include/iterative_solvers.hpp"

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
  int N = 200;
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



  //resvec + time
  A = create_A_N(N);
  b.resize(N);
  b.fill(1.0);
  vector<double> resvec;
  file.open("resvec.json");

  file << "{\n\t";

  //Jacobi
  file << "\"Jacobi\" : {\n\t\t\"N\" : " << N << ",\n\t\t\"resvec\" : [";
  Jacobi_.solve();
  resvec = Jacobi_.get_resvec();
  file << resvec[0];
  for(unsigned int i=1;i<resvec.size();i++){
    file << ", " << resvec[i];
  }
  file << "],\n\t\t\"sizes\" : [";
  for(auto ptr = sizes.begin();ptr!=sizes.end();ptr++){
    if(ptr!=sizes.begin()){
      file << ", ";
    }
    file << *ptr;
  }
  file << "],\n\t\t\"duration\" : [";
  for(auto ptr = sizes.begin();ptr!=sizes.end();ptr++){
    A = create_A_N((*ptr));
    b.resize((*ptr));
    b.fill(1.0);
    Jacobi_.solve();
    if(ptr!=sizes.begin()){
      file << ", ";
    }
    file << Jacobi_.get_duration();
  }
  file << "]\n\t},\n";

  //GaussSeidel
  file << "\n\t\"GaussSeidel\" : {\n\t\t\"N\" : " << N << ",\n\t\t\"resvec\" : [";
  GaussSeidel_.solve();
  resvec = GaussSeidel_.get_resvec();
  file << resvec[0];
  for(unsigned int i=1;i<resvec.size();i++){
    file << ", " << resvec[i];
  }
  file << "],\n\t\t\"sizes\" : [";
  for(auto ptr = sizes.begin();ptr!=sizes.end();ptr++){
    if(ptr!=sizes.begin()){
      file << ", ";
    }
    file << *ptr;
  }
  file << "],\n\t\t\"duration\" : [";
  for(auto ptr = sizes.begin();ptr!=sizes.end();ptr++){
    A = create_A_N((*ptr));
    b.resize((*ptr));
    b.fill(1.0);
    GaussSeidel_.solve();
    if(ptr!=sizes.begin()){
      file << ", ";
    }
    file << GaussSeidel_.get_duration();
  }
  file << "]\n\t},\n";

  //Relax
  file << "\n\t\"Relax\" : {\n\t\t\"N\" : " << N << ",\n\t\t\"resvec\" : [";
  rho = cos(M_PI/(N+1));
  omg_star = 2.0/(1.0+sqrt( 1.0 - rho*rho ));
  Relax_.set_omega(omg_star);
  Relax_.solve();
  resvec = Relax_.get_resvec();
  file << resvec[0];
  for(unsigned int i=1;i<resvec.size();i++){
    file << ", " << resvec[i];
  }
  file << "],\n\t\t\"sizes\" : [";
  for(auto ptr = sizes.begin();ptr!=sizes.end();ptr++){
    if(ptr!=sizes.begin()){
      file << ", ";
    }
    file << *ptr;
  }
  file << "],\n\t\t\"duration\" : [";
  for(auto ptr = sizes.begin();ptr!=sizes.end();ptr++){
    A = create_A_N((*ptr));
    b.resize((*ptr));
    b.fill(1.0);
    rho = cos(M_PI/((*ptr)+1));
    omg_star = 2.0/(1.0+sqrt( 1.0 - rho*rho ));
    Relax_.set_omega(omg_star);
    Relax_.solve();
    if(ptr!=sizes.begin()){
      file << ", ";
    }
    file << Relax_.get_duration();
  }
  file << "]\n\t}\n";


  file << "}";
  file.close();



  return 0;
}
