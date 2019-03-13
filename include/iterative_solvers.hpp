#ifndef ITERATIVE_SOLVERS
#define ITERATIVE_SOLVERS

#include <chrono>

#include "direct_solvers.hpp"

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
  std::chrono::high_resolution_clock::time_point begin, end;

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
  const Vect& get_x();
  int get_niter();
  const std::vector<double>& get_resvec();
  const Vect& get_r();
  long int get_duration();

  void solve();
  virtual void check();
  virtual void update_solution() = 0;
  virtual void update_resvec();

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

class GaussSeidel final : public IterSolver{

  static GaussSeidel obj;

  GaussSeidel();
  GaussSeidel(const GaussSeidel&) = delete;
  GaussSeidel(GaussSeidel&&) = delete;
  GaussSeidel& operator=(const GaussSeidel&) = delete;
  GaussSeidel& operator=(GaussSeidel&&) = delete;

public:

  void check() override;
  void update_solution() override;
  static GaussSeidel& getobj();

};

class Relax final : public IterSolver{

  static Relax obj;
  double omega_;

  Relax();
  Relax(const Relax&) = delete;
  Relax(Relax&&) = delete;
  Relax& operator=(const Relax&) = delete;
  Relax& operator=(Relax&&) = delete;

public:

  void check() override;
  void update_solution() override;
  static Relax& getobj();
  void set_omega(double);

};


extern Jacobi& Jacobi_;
extern GaussSeidel& GaussSeidel_;
extern Relax& Relax_;

#endif
