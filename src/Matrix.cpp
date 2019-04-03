#include "../include/Matrix.hpp"
#include "../include/direct_solvers.hpp"
#include <cmath>
#include <cassert>
#include <typeinfo>
#include <algorithm>


using namespace std;


Matrix::Matrix(int N) : N_(N){}

int Matrix::size() const{
  return N_;
}


Vect operator*(const Matrix& A, const Vect& v){
  if(typeid(A)==typeid(DenseMatrix)){
    return (dynamic_cast<const DenseMatrix&>(A))*v;
  }
  return (dynamic_cast<const COO&>(A))*v;;
}





DenseMatrix::DenseMatrix() : Matrix(0), coef_(0){}

DenseMatrix::DenseMatrix(int N) : Matrix(N), coef_(N*N){}


void DenseMatrix::display() const{
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

void DenseMatrix::display_LU() const{

  if(state!=LU){
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

void DenseMatrix::display_Cholesky() const{

  if(state!=CHOLESKY){
    cout << "The Cholesky decomposition has not been computed\n";
    return;
  }
  cout << scientific << setprecision(2) << "N_ : " << N_ << "\n";
  cout << "A = L(L^T) with\n";
  //L
  for(int i=0;i<N_;i++){
    if(N_/2==i){
      cout << "L = |";
    }else{
      cout << "    |";
    }
    for(int j=0;j<=i;j++){
      cout << setw(10) << coef_[i*N_+j];
    }
    for(int j=i+1;j<N_;j++){
      cout << setw(10) << 0.0;
    }
    cout << "|\n";
  }
  cout << "\n\n";
}

void DenseMatrix::diag(double alpha, int d){
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

void DenseMatrix::decomp_LU(){

  int i, j, k, kmax;
  double aux;

  if(state!=NORMAL){
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

  state = LU;
}

void DenseMatrix::cholesky(){
  int i, j, k;

  if(state!=NORMAL){
    return;
  }

  for(k=0;k<N_;k++){
    coef_[k*(N_+1)] = sqrt(coef_[k*(N_+1)]);
    for(i=1;i<N_-k;i++){
      coef_[(k+i)*N_+k] /= coef_[k*(N_+1)];
      coef_[k*N_+k+i] = coef_[(k+i)*N_+k];
      for(j=1;j<i;j++){
	coef_[(k+i)*N_+k+j] -= coef_[(k+i)*N_+k]*coef_[(k+j)*N_+k];
	coef_[(k+j)*N_+k+i] = coef_[(k+i)*N_+k+j];
      }
      coef_[(k+i)*(N_+1)] -= coef_[(k+i)*N_+k]*coef_[(k+i)*N_+k];
    }
  }

  state = CHOLESKY;
}

Vect DenseMatrix::solve_via_LU(const Vect& b){

  assert(state==LU);

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

Vect DenseMatrix::solve_via_Cholesky(const Vect& b){

  assert(state==CHOLESKY);

  //we solve LL^T x = b
  Vect c(solve_triang_inf(*this, b));
  return solve_triang_sup(*this, c);
}

double& DenseMatrix::operator()(int i, int j){
  return coef_[i*N_+j];
}

double DenseMatrix::operator()(int i, int j) const{
  return coef_[i*N_+j];
}


ostream& operator<<(ostream& os, const DenseMatrix& m){
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

DenseMatrix operator*(double alpha, const DenseMatrix& m){
  DenseMatrix res(m.N_);
  for(unsigned int i=0;i<m.coef_.size();i++){
    res.coef_[i] = alpha*m.coef_[i];
  }
  return res;
}

DenseMatrix operator+(const DenseMatrix& A, const DenseMatrix& B){

  if(A.N_!=B.N_){
    throw(logic_error("DenseMatrix operator+(const DenseMatrix& A, const DenseMatrix& B) : matrices must have the same size\n"));
  }

  unsigned int i;
  DenseMatrix res(A.N_);

  for(i=0;i<A.coef_.size();i++){
    res.coef_[i] = A.coef_[i] + B.coef_[i];
  }

  return res;
}

DenseMatrix operator-(const DenseMatrix& A, const DenseMatrix& B){

  if(A.N_!=B.N_){
    throw(logic_error("DenseMatrix operator-(const DenseMatrix& A, const DenseMatrix& B) : matrices must have the same size\n"));
  }

  unsigned int i;
  DenseMatrix res(A.N_);

  for(i=0;i<A.coef_.size();i++){
    res.coef_[i] = A.coef_[i] - B.coef_[i];
  }

  return res;
}

DenseMatrix operator*(const DenseMatrix& A, const DenseMatrix& B){

  if(A.N_!=B.N_){
    throw(logic_error("DenseMatrix operator*(const DenseMatrix& A, const DenseMatrix& B) : matrices must have the same size\n"));
  }

  int i, j, k;
  DenseMatrix res(A.N_);

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

Vect operator*(const DenseMatrix& A, const Vect& v){

  if(A.N_!=v.size()){
    throw(logic_error("DenseMatrix operator*(const DenseMatrix& A, const Vect& v) : the matrice A and the vector v must have the same size\n"));
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







COO::COO(const DenseMatrix& A) : Matrix(A.size()), data(0){
  data.reserve(A.size()*A.size());
  int i, j;
  tuple3<int, int, double> aux;
  for(i = 0;i<A.size();i++){
    aux.first = i;
    for(j = 0;j<A.size();j++){
      if(A(i, j)!=0.0){
	aux.second = j;
	aux.third = A(i, j);
	data.push_back(aux);
      }
    }
  }
  data.shrink_to_fit();
}

COO COO::Laplacian(int n){
  COO Lapl(n);
  int i, j;
  tuple3<int, int, double> aux;

  int length;
  switch(n){
  case 1:
    length = 1;
    aux.first = 0;
    aux.second = 0;
    aux.third = 2.0;
    Lapl.data.push_back(aux);
    return Lapl;
  case 2:
    length = 4;
    break;
  default:
    length = (n-2)*3+4;
  }

  Lapl.data.reserve(length);

  aux.first = 0;
  aux.second = 0;
  aux.third = 2.0;
  Lapl.data.push_back(aux);
  aux.second = 1;
  aux.third = -1.0;
  Lapl.data.push_back(aux);
  for(i = 1;i<n-1;i++){
    aux.first = i;
    for(j = -1;j<=1;j++){
      aux.second = i+j;
      aux.third = (j==0)?2.0:-1.0;
      Lapl.data.push_back(aux);
    }
  }
  aux.first = n-1;
  aux.second = n-2;
  aux.third = -1.0;
  Lapl.data.push_back(aux);
  aux.second = n-1;
  aux.third = 2.0;
  Lapl.data.push_back(aux);

  return Lapl;
}

double& COO::operator()(int i, int j){
  for(tuple3<int, int, double>& d : data){
    if(d.first==i && d.second==j){
      return d.third;
    }
  }
  tuple3<int, int, double> aux;
  aux.first = i;
  aux.second = j;
  data.push_back(aux);
  return (data.end()-1)->third;
}

double COO::operator()(int i, int j) const{
  for(const tuple3<int, int, double>& d : data){
    if(d.first==i && d.second==j){
      return d.third;
    }
  }
  return 0.0;
}

void COO::display() const{
  cout << "N_ : " << N_ << "\ndata :\n";
  for(const tuple3<int, int, double>& t : data){
    cout << "\t( " << t.first << ", " << t.second << ", " << setw(10) << t.third << " )\n";
  }
}

int COO::get_nnz() const{
  return data.size();
}

Vect operator*(const COO& A, const Vect& v){
  Vect res(v.size());
  for(int i=0;(vector<tuple3<int, int, double> >::size_type)i<A.data.size();i++){
    res(A.data[i].first) += A.data[i].third * v(A.data[i].second);
  }
  return res;
}



CSR::CSR(const COO& A) : Matrix(A.size()), row(0), col(0), val(0){
  //A is supposed without null lign
  row.reserve(A.size());
  col.reserve(A.get_nnz());
  val.reserve(A.get_nnz());

  sort(A.data.begin(), A.data.end());

  row.push_back(0);
  col.push_back(A.data[0].second);
  val.push_back(A.data[0].third);
  for(int i=1;i<(int)(A.data.size());i++){
    if(A.data[i].second>A.data[i-1].second){
      row.push_back(col.size());
    }
    col.push_back(A.data[i].second);
    val.push_back(A.data[i].third);
  }
}


double& operator()(int, int){

}

double operator()(int, int) const{

}

Vect operator*(const CSR&, const Vect&){

}
