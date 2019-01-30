#ifndef ITERATIVE_SOLVERS
#define ITERATIVE_SOLVERS

#include "Matrix.hpp"

class IterSolver{

  //input
  Matrix A_;
  Vect b_;
  double tol_;
  int n_max_;

  //output
  Vect x;
  int niter_;
  std::vector<double> resvec_;


public:






};

class Jacobi : public IterSolver{

public:

};

#endif
