#include "../include/direct_solvers.hpp"

using namespace std;



Vect solve_triang_inf(const Matrix& A, const Vect& b){

  Vect res(b);
  int i, j;

  for(i = 0;i<b.size();i++){
    for(j = 0;j<i;j++){
      res(i) -= A(i, j)*res(j);
    }
    res(i) /= A(i, i);
  }

  return res;
}

Vect solve_triang_sup(const Matrix& A, const Vect& b){

  Vect res(b);
  int i, j;

  for(i = b.size()-1;i>=0;i--){
    for(j = i+1;j<b.size();j++){
      res(i) -= A(i, j)*res(j);
    }
    res(i) /= A(i, i);
  }

  return res;

}

Vect solve_triang_inf_id(const Matrix& A, const Vect& b){
  Vect res(b);
  int i, j;

  for(i = 0;i<b.size();i++){
    for(j = 0;j<i;j++){
      res(i) -= A(i, j)*res(j);
    }
  }

  return res;
}

Vect solve_triang_sup_id(const Matrix& A, const Vect& b){

  Vect res(b);
  int i, j;

  for(i = b.size()-1;i>=0;i--){
    for(j = i+1;j<b.size();j++){
      res(i) -= A(i, j)*res(j);
    }
  }

  return res;

}
