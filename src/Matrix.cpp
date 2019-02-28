#include "../include/Matrix.hpp"
#include "../include/direct_solvers.hpp"
#include <cmath>

using namespace std;

Matrix::Matrix() : N_(3), coef_(9){}

Matrix::Matrix(int N) : N_(N), coef_(N*N){}


int Matrix::size() const{
  return N_;
}

void Matrix::display() const{
  cout << scientific << setprecision(2) << "N_ : " << N_ << "\ncoef_ :\n";
  for(unsigned int i=0;i<coef_.size();i++){
    if(i%((unsigned int)N_)==0u){
      cout << "|";
    }
    cout << " " << setw(10) << coef_[i];
    if((i%((unsigned int)N_))+1u==(unsigned int)N_){
      cout << " |\n";
    }
  }
  cout << "\n\n";
}

void Matrix::display_LU() const{

  if(!LU_){
    cout << "The LU decomposition has not been computed\n";
    return;
  }
  cout << scientific << setprecision(2) << "N_ : " << N_ << "\n";
  cout << "PA = LU with\n";
  //L
  for(int i=0;i<N_;i++){
    if(N_/2==i){
      cout << "L = |";
    }else{
      cout << "    |";
    }
    for(int j=0;j<i;j++){
      cout << setw(10) << coef_[i*N_+j];
    }
    cout << setw(10) << 1.0;
    for(int j=i+1;j<N_;j++){
      cout << setw(10) << 0.0;
    }
    cout << "|\n";
  }
  cout << "\n\n";
  //U
  for(int i=0;i<N_;i++){
    if(N_/2==i){
      cout << "U = |";
    }else{
      cout << "    |";
    }
    for(int j=0;j<i;j++){
      cout << setw(10) << 0.0;
    }
    for(int j=i;j<N_;j++){
      cout << setw(10) << coef_[i*N_+j];
    }
    cout << "|\n";
  }

  //perm
  cout << "P is given by those operations (in that order):\n";
  for(int i=0;i<N_-1;i++){
    cout << "\t" << i << " <-> " << perm_[i] << "\n";
  }
}

void Matrix::diag(double alpha, int d){
  if(d>=0){
    for(int i=d;i<N_;i++){
      coef_[i*N_+i-d] = alpha;
    }
  }else{
    for(int i=0;i<N_+d;i++){
      coef_[i*N_+i-d] = alpha;
    }
  }
}

void Matrix::decomp_LU(){

  int i, j, k, kmax;
  double aux;

  if(LU_){
    return;
  }

  perm_.reserve(N_-1);

  for(k=0;k<=N_-2;k++){

    //we select the maximum of the entries below the entry (k, k) (included)
    //then we swap two rows
    aux = abs(coef_[k*(N_+1)]);
    kmax = k;
    for(i=1;i<N_-k;i++){
      if(abs(coef_[(k+i)*N_+k])>aux){
	aux = abs(coef_[(k+i)*N_+k]);
	kmax = k+i;
      }
    }
    perm_[k] = kmax;
    for(i=k;i<N_;i++){
      aux = coef_[k*N_+i];
      coef_[k*N_+i] = coef_[kmax*N_+i];
      coef_[kmax*N_+i] = aux;
    }

    for(i=1;i<N_-k;i++){
      coef_[(k+i)*N_+k] /= coef_[k*(N_+1)];
      for(j=1;j<N_-k;j++){
	coef_[(k+i)*N_+k+j] -= coef_[(k+i)*N_+k]*coef_[k*N_+j+k];
      }
    }
  }

  LU_ = true;
}

Vect Matrix::solve_via_LU(const Vect& b){
  if(!LU_){
    this->decomp_LU();
  }
  //we solve LUx = Pb
  Vect c(b);
  double aux;
  for(int i=0;i<N_-1;i++){
    aux = c(i);
    c(i) = c(perm_[i]);
    c(perm_[i]) = aux;
  }
  //c is ready and equals Pb
  c = solve_triang_inf_id(*this, c);
  return solve_triang_sup(*this, c);
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
