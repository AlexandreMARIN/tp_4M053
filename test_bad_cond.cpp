#include "./include/iterative_solvers.hpp"
#include <cmath>
#include <fstream>

using namespace std;


Matrix create_bad_matrix(int n, double alpha){
  Matrix A(n);
  int i, j;
  for(i=0;i<n;i++){
    for(j=0;j<i;j++){
      A(i, j) = 1.0/(1.0 + pow(abs(i-j), alpha));
      A(j, i) = A(i, j);
    }
    A(i, i) = 1.0;
  }

  return A;
}


int main(){

  int n;

  cout << "Give an integer : ";
  cin >> n;

  if(cin.bad()){
    cout << "\nn = 2\n";
    n = 2;
  }

  ofstream file("bad_matrices.json");
  Matrix A(0);
  Vect b(n);
  int n_max = 1000;
  double tol = 1e-3;
  b.fill(1.0);
  vector<double> resvec;
  vector<double> alpha{2.0, 1.5, 1.2, 1.1, 1.05};
  Relax_.set_omega(1.5);
  Relax_.set_tol(tol);
  Relax_.set_n_max(n_max);
  Relax_.set_A(&A);
  Relax_.set_b(&b);

  file << "{\n\t\"N\" : " << n << ",\n\t\"alpha\" : [";
  file << alpha[0];
  for(auto ptr = alpha.begin()+1;ptr!=alpha.end();ptr++){
    file << ", " << *ptr;
  }

  file << "],\n";
  file << "\t\"Relax\" : {\n\t\t\"resvec\" : [";

  for(auto a = alpha.begin();a!=alpha.end();a++){
    A = create_bad_matrix(n, *a);
    Relax_.solve();
    resvec = Relax_.get_resvec();
    file << "[" << resvec[0];
    for(auto r = resvec.begin()+1;r!=resvec.end();r++){
      file << ", " << *r;
    }
    file << "]";
    if(a!=alpha.end()-1){
      file << ",";
    }
  }

  file << "]\n\t}\n}";

  return 0;
}
