#include "../include/conj_grad.hpp"

using namespace std;

ConjGrad ConjGrad::obj;

ConjGrad::ConjGrad() : IterSolver(){}

void ConjGrad::check(){
  p_ = *b_;//p_0 = r_0
  //n_max_ = b_->size();
}

volatile void ConjGrad::update_solution(){
  aux = (r_, r_);
  alpha_ = aux/(*A_*p_, p_);
  x_ += alpha_*p_;
}

volatile void ConjGrad::update_resvec(){
  r_ -= alpha_*((*A_)*p_);
  aux = (r_, r_)/aux;
  p_ = r_ + aux*p_;
}

ConjGrad& ConjGrad::getobj(){
  return obj;
}

ConjGrad& ConjGrad_ = ConjGrad::getobj();
