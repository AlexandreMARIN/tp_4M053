#ifndef OPT_STEP_GRAD_HPP
#define OPT_STEP_GRAD_HPP

#include "iterative_solvers.hpp"

class OptStepGrad final: public IterSolver{


  double alpha_;
  Vect aux;

  OptStepGrad();
  OptStepGrad(const OptStepGrad&) = delete;
  OptStepGrad(OptStepGrad&&) = delete;
  OptStepGrad& operator=(const OptStepGrad&) = delete;
  OptStepGrad& operator=(OptStepGrad&&) = delete;

public:

  static OptStepGrad obj;

  void check() override;
  volatile void update_solution() override;
  volatile void update_resvec() override;
  static OptStepGrad& getobj();

};

extern OptStepGrad& OptStepGrad_;

#endif
