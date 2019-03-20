#include <fstream>
#include <cmath>


#include "include/conj_grad.hpp"
#include "include/opt_step_grad.hpp"

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
  constexpr int measures = 3;

  vector<int> sizes(measures);
  sizes[0] = 20;
  for(int i=1;i<measures;i++){
    sizes[i] = sizes[i-1]*2;
  }

  double tol = 1e-2;
  int n_max = 1000;
  int N = 50;
  Matrix A(0);
  Vect b(0);

  double omg_star, rho;

  cout << "Give an integer for N : ";
  cin >> N;

  if(cin.bad()){
    cout << "\nN = 50\n";
    N = 50;
  }


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

  OptStepGrad_.set_tol(tol);
  OptStepGrad_.set_n_max(n_max);
  OptStepGrad_.set_A(&A);
  OptStepGrad_.set_b(&b);

  ConjGrad_.set_tol(tol);
  ConjGrad_.set_n_max(n_max);
  ConjGrad_.set_A(&A);
  ConjGrad_.set_b(&b);

  //resvec + time
  A = create_A_N(N);
  b.resize(N);
  b.fill(1.0);
  vector<double> resvec;
  file.open("resvec2.json");

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

  //OptStepGrad
  file << "\n\t\"OptStepGrad\" : {\n\t\t\"N\" : " << N << ",\n\t\t\"resvec\" : [";
  OptStepGrad_.solve();
  resvec = OptStepGrad_.get_resvec();
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
    OptStepGrad_.solve();
    if(ptr!=sizes.begin()){
      file << ", ";
    }
    file << OptStepGrad_.get_duration();
  }
  file << "]\n\t},\n";

  //ConjGrad
  file << "\n\t\"ConjGrad\" : {\n\t\t\"N\" : " << N << ",\n\t\t\"resvec\" : [";
  ConjGrad_.solve();
  resvec = ConjGrad_.get_resvec();
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
    ConjGrad_.solve();
    if(ptr!=sizes.begin()){
      file << ", ";
    }
    file << ConjGrad_.get_duration();
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
