#ifndef ITERATIVE_SOLVERS
#define ITERATIVE_SOLVERS

#include "Matrix.hpp"

class IterSolver{
protected:
  //input
  Matrix* A_;
  Vect* b_;
  double tol_;
  int n_max_;

  //output
  Vect x_;
  int niter_;
  std::vector<double> resvec_;

  Vect r_;


  IterSolver();
  IterSolver(const IterSolver&) = delete;
  IterSolver(IterSolver&&) = delete;
  IterSolver& operator=(const IterSolver&) = delete;
  IterSolver& operator=(IterSolver&&) = delete;

public:

  void set_A(Matrix*);
  void set_b(Vect*);
  void set_tol(double);
  void set_n_max(int);
  Vect& get_x();
  void solve();
  virtual void check();
  virtual void update_solution() = 0;

};

class Jacobi final : public IterSolver{

  static Jacobi obj;

  Jacobi();
  Jacobi(const Jacobi&) = delete;
  Jacobi(Jacobi&&) = delete;
  Jacobi& operator=(const Jacobi&) = delete;
  Jacobi& operator=(Jacobi&&) = delete;

public:

  void check() override;
  void update_solution() override;
  static Jacobi& getobj();

};

extern Jacobi& Jacobi_;

#endif
