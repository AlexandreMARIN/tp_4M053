#ifndef CONJ_GRAD_HPP
#define CONJ_GRAD_HPP

#include "iterative_solvers.hpp"

class ConjGrad final: public IterSolver{


  double alpha_;
  double aux;
  Vect p_;

  ConjGrad();
  ConjGrad(const ConjGrad&) = delete;
  ConjGrad(ConjGrad&&) = delete;
  ConjGrad& operator=(const ConjGrad&) = delete;
  ConjGrad& operator=(ConjGrad&&) = delete;

public:

  static ConjGrad obj;

  void check() override;
  volatile void update_solution() override;
  volatile void update_resvec() override;
  static ConjGrad& getobj();

};

extern ConjGrad& ConjGrad_;

#endif
