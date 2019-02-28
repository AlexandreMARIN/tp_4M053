#include <fstream>
#include <cmath>
#include <map>
#include <string>


#include "include/iterative_solvers.hpp"

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
  map<string, vector<vector<int> > > data;//use pairs or somethong like that
  double tol = 0.1;
  int n_max = 100000;

  Jacobi_.set_tol(tol);
  Jacobi_.set_n_max(n_max);

  GaussSeidel_.set_tol(tol);
  GaussSeidel_.set_n_max(n_max);

  Relax_.set_tol(tol);
  Relax_.set_n_max(n_max);



  ofstream file("iter.json");
  Matrix A(0);
  Vect b(0);
  vector<string> methods{"Jacobi", "GaussSeidel", "Relax"};


  Jacobi_.set_A(&A);
  Jacobi_.set_b(&b);

  GaussSeidel_.set_A(&A);
  GaussSeidel_.set_b(&b);

  Relax_.set_A(&A);
  Relax_.set_b(&b);

  for(string& meth : methods){
    data[meth] = vector<vector<int> >(3);
    for(int i=0;i<3;i++){
      data.at(meth)[i].reserve(sizes.size());
    }
  }

  for(int s : sizes){
    A = create_A_N(s);
    b.resize(s);
    b.fill(1.0);

    //Jacobi
    Jacobi_.solve();
    data.at("Jacobi")[0].push_back(Jacobi_.get_niter());

    //"analytical assessment"
    data.at("Jacobi")[1].push_back( ceil(log(tol)/log(cos(M_PI/(s+1)))) );

    //assessment no 3
    data.at("Jacobi")[2].push_back( -2.0*log(tol)*(s+1)*(s+1)/(M_PI*M_PI) );

    //GaussSeidel
    GaussSeidel_.solve();
    data.at("GaussSeidel")[0].push_back( GaussSeidel_.get_niter() );

    //"analytical assessment"
    data.at("GaussSeidel")[1].push_back( ceil(log(tol)/(2.0*log(cos(M_PI/(s+1))))) );

    //assessment no 3
    data.at("GaussSeidel")[2].push_back( ceil(-log(tol)*(s+1)*(s+1)/(M_PI*M_PI)) );

    //Relax
    double rho = cos(M_PI/(s+1));
    double omg_star = 2.0/(1.0+sqrt( 1.0 - rho*rho ));
    Relax_.set_omega(omg_star);

    Relax_.solve();
    data.at("Relax")[0].push_back( Relax_.get_niter() );

    //"analytical assessment"
    data.at("Relax")[1].push_back( ceil(log(tol)/log(omg_star - 1.0)) );

    //assessment no 3
    data.at("Relax")[2].push_back( ceil(-log(tol)*(s+1)/(M_PI*M_PI)) );
  }


  //we write the file here
  //we use Json
  vector<string> legend{"niter", "approx1", "approx2"};
  file << "{\n";
  //we write sizes, then measures are written
  file << "\t\"sizes\" : [";
  for(auto si=sizes.begin();si!=sizes.end();si++){
    file << *si;
    if(si!=sizes.end()-1){
      file << ", ";
    }
  }
  file << "],\n\n";

  for(auto meth=methods.begin();meth!=methods.end();meth++){
    file << "\t\"" << *meth << "\" : {\n\t\t";
    for(int i=0;i<3;i++){
      file << "\"" << legend[i] << "\" : [";
      for(auto ptr=data.at(*meth)[i].begin();ptr!=data.at(*meth)[i].end();ptr++){
	file << *ptr;
	if(ptr!=data.at(*meth)[i].end()-1){
	  file << ", ";
	}
      }
      file << "]";
      if(i<2){
	file << ",\n\t\t";
      }else{
	file << "\n\t}";
      }
    }
    if(meth!=methods.end()-1){
      file << ",";
    }
    file << "\n";
  }
  file << "\n}";
  file.close();


  return 0;
}
