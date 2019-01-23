#include "../include/Matrix.hpp"

using namespace std;

Matrix::Matrix() : N_(3), coef_(9){}

Matrix::Matrix(int N) : N_(N), coef_(N*N){}


int Matrix::size() const{
  return N_;
}

void Matrix::display() const{
  cout << "N_ : " << N_ << "\ncoef_ :\n";
  for(unsigned int i=0;i<coef_.size();i++){
    if(i%((unsigned int)N_)==0u){
      cout << "|";
    }
    cout << " " << coef_[i];
    if((i%((unsigned int)N_))+1u==(unsigned int)N_){
      cout << " |\n";
    }
  }
  cout << "\n\n";
}

double& Matrix::operator()(int i, int j){
  return coef_[i*N_+j];
}

double Matrix::operator()(int i, int j) const{
  return coef_[i*N_+j];
}


ostream& operator<<(ostream& os, const Matrix& m){
  os << "N_ : " << m.N_ << "\ncoef_ :\n";
  for(unsigned int i=0;i<m.coef_.size();i++){
    if(i%((unsigned int)m.N_)==0u){
      os << "|";
    }
    os << " " << m.coef_[i];
    if((i%((unsigned int)m.N_))+1u==(unsigned int)m.N_){
      os << " |\n";
    }
  }
  os << "\n\n";
  return os;
}

Matrix operator*(double alpha, const Matrix& m){
  Matrix res(m.N_);
  for(unsigned int i=0;i<m.coef_.size();i++){
    res.coef_[i] = alpha*m.coef_[i];
  }
  return res;
}

Matrix operator+(const Matrix& A, const Matrix& B){

  if(A.N_!=B.N_){
    throw(logic_error("Matrix operator+(const Matrix& A, const Matrix& B) : matrices must have the same size\n"));
  }

  unsigned int i;
  Matrix res(A.N_);

  for(i=0;i<A.coef_.size();i++){
    res.coef_[i] = A.coef_[i] + B.coef_[i];
  }

  return res;
}

Matrix operator-(const Matrix& A, const Matrix& B){

  if(A.N_!=B.N_){
    throw(logic_error("Matrix operator-(const Matrix& A, const Matrix& B) : matrices must have the same size\n"));
  }

  unsigned int i;
  Matrix res(A.N_);

  for(i=0;i<A.coef_.size();i++){
    res.coef_[i] = A.coef_[i] - B.coef_[i];
  }

  return res;
}

Matrix operator*(const Matrix& A, const Matrix& B){

  if(A.N_!=B.N_){
    throw(logic_error("Matrix operator*(const Matrix& A, const Matrix& B) : matrices must have the same size\n"));
  }

  int i, j, k;
  Matrix res(A.N_);

  for(i=0;i<A.N_;i++){
    for(j=0;j<A.N_;j++){
      res.coef_[i*A.N_+j] = 0.0;
      for(k=0;k<A.N_;k++){
	res.coef_[i*A.N_+j] += A.coef_[i*A.N_+k] * B.coef_[k*A.N_+j];
      }
    }
  }

  return res;
}

Vect operator*(const Matrix& A, const Vect& v){

  if(A.N_!=v.size()){
    throw(logic_error("Matrix operator*(const Matrix& A, const Vect& v) : the matrice A and the vector v must have the same size\n"));
  }

  int i, j;
  Vect res(v.size());

  for(i=0;i<A.N_;i++){
    res(i) = 0.0;
    for(j=0;j<A.N_;j++){
      res(i) += A.coef_[i*A.N_+j] * v(j);
    }
  }

  return res;
}
