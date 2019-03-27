#include <limits>
#include "../include/iterative_solvers.hpp"

using namespace std;
using namespace std::chrono;

IterSolver::IterSolver() : A_(nullptr), b_(nullptr), tol_(1e-6), n_max_(100){}

void IterSolver::set_A(Matrix* A){
  A_ = A;
}

void IterSolver::set_b(Vect* b){
  b_ = b;
}

void IterSolver::set_tol(double tol){
  if(tol>0.){
    tol_ = tol;
  }
}

void IterSolver::set_n_max(int n_max){
  if(n_max>0){
    n_max_ = n_max;
  }
}

const Vect& IterSolver::get_x(){
  return x_;
}

int IterSolver::get_niter(){
  return niter_;
}

const vector<double>& IterSolver::get_resvec(){
  return resvec_;
}

const Vect& IterSolver::get_r(){
  return r_;
}

long int IterSolver::get_duration(){
  milliseconds dur = duration_cast<milliseconds>(end - begin);
  return dur.count();
}

void IterSolver::solve(){
  IterSolver::check();
  this->check();
  double rnorm, bnorm=b_->norm_infty(), limit=bnorm*tol_;

  resvec_.resize(0);
  resvec_.reserve(n_max_+1);
  niter_ = 0;
  x_.resize(b_->size());
  x_.set_to_zero();
  r_ = *b_;

  begin = high_resolution_clock::now();
  while( (rnorm=r_.norm_infty())>limit && niter_<n_max_ ){
    update_solution();
    this->update_resvec();
    resvec_.push_back((bnorm!=0.0)?rnorm/bnorm:numeric_limits<double>::quiet_NaN());
    niter_++;
  }
  end = high_resolution_clock::now();
  resvec_.push_back((bnorm!=0.0)?rnorm/bnorm:numeric_limits<double>::quiet_NaN());

}

void IterSolver::check(){

  if(!A_ || !b_){
    throw(invalid_argument("IterSolver : members A_ and b_ are set to null\n"));
  }

  /*  if(A_->size()!=b_->size()){
    throw(invalid_argument("IterSolver : matrix A_ and vector b_ must have the same size\n"));
    }*/

}

volatile void IterSolver::update_resvec(){
    r_ = *b_ - (*A_)*x_;
}


Jacobi Jacobi::obj;

Jacobi::Jacobi() : IterSolver(){}

void Jacobi::check(){

  //A_'s diagonal must not contain any zero
  /*  for(int i=0;i<A_->size();i++){
    if((*A_)(i, i)==0.0){
      throw(invalid_argument("Jacobi : A_'s diagonal must not contain any zero\n"));
    }
    }*/

}

volatile void Jacobi::update_solution(){

  for(int i=0;i<x_.size();i++){
    x_(i) += r_(i)/(*A_)(i, i);
  }

}

Jacobi& Jacobi::getobj(){
  return obj;
}



GaussSeidel GaussSeidel::obj;

GaussSeidel::GaussSeidel() : IterSolver(){}

void GaussSeidel::check(){

  //A_'s diagonal must not contain any zero
  /*  for(int i=0;i<A_->size();i++){
    if((*A_)(i, i)==0.0){
      throw(invalid_argument("GaussSeidel : A_'s diagonal must not contain any zero\n"));
    }
    }*/

}

volatile void GaussSeidel::update_solution(){

  Vect y = solve_triang_sup(*A_, r_);
  x_ += y;

}

GaussSeidel& GaussSeidel::getobj(){
  return obj;
}


Relax Relax::obj;

Relax::Relax() : IterSolver(), omega_(0.0){}

void Relax::check(){

  if(omega_<=0. || omega_>=2.0){
    throw(invalid_argument("Relax : omega_ must be in ]0, 2[\n"));
  }

  //A_'s diagonal must not contain any zero
  /*  for(int i=0;i<A_->size();i++){
    if((*A_)(i, i)==0.0){
      throw(invalid_argument("Relax : A_'s diagonal must not contain any zero\n"));
    }

    }*/

}

volatile void Relax::update_solution(){

  Vect y(r_);
  int i, j;

  for(i = r_.size()-1;i>=0;i--){
    for(j = i+1;j<r_.size();j++){
      y(i) -= (*A_)(i, j)*y(j);
    }
    y(i) /= (*A_)(i, i);
    y(i) *= omega_;
  }

  x_ += y;

}

Relax& Relax::getobj(){
  return obj;
}

void Relax::set_omega(double o){
  omega_ = o;
}

Jacobi& Jacobi_ = Jacobi::getobj();
GaussSeidel& GaussSeidel_ = GaussSeidel::getobj();
Relax& Relax_ = Relax::getobj();
