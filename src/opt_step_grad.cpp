#include "../include/opt_step_grad.hpp"

using namespace std;

OptStepGrad OptStepGrad::obj;

OptStepGrad::OptStepGrad() : IterSolver(){}

void OptStepGrad::check(){}

void OptStepGrad::update_solution(){
  double normA = (r_,((*A_)*r_));
  if(normA==0.0){
    n_max_ = niter_;
    return;
  }
  alpha_ = (r_,r_)/normA;
  aux = alpha_*r_;
  x_ += aux;
}

void OptStepGrad::update_resvec(){
  r_ -= (*A_)*aux;
}

OptStepGrad& OptStepGrad::getobj(){
  return obj;
}

OptStepGrad& OptStepGrad_ = OptStepGrad::getobj();
