#ifndef MATRIX_HPP
#define MATRIX_HPP

#include <vector>
#include <iostream>


#include "Vect.hpp"


class Matrix{
public:
  Matrix() = default;
  Matrix(const Matrix&) = default;
  Matrix(Matrix&&) = default;

  Matrix& operator=(const Matrix&) = default;
  Matrix& operator=(Matrix&&) = default;

  virtual double& operator()(int, int) = 0;
  virtual double operator()(int, int) const = 0;

};

Vect operator*(const Matrix&, const Vect&);


class DenseMatrix : public Matrix{
private:
  enum state_ {NORMAL, LU, CHOLESKY} state = NORMAL;

  int N_;//number of lines and number of columns
  std::vector<double> coef_;//entries
  bool LU_ = false;
  bool cholesky_ = false;
  std::vector<int> perm_;
public:
  DenseMatrix();
  DenseMatrix(int);
  DenseMatrix(const DenseMatrix&) = default;
  DenseMatrix(DenseMatrix&&) = default;
  ~DenseMatrix() = default;

  int size() const;
  void display() const;
  void display_LU() const;
  void display_Cholesky() const;

  void diag(double, int=0);
  void decomp_LU();
  void cholesky();
  Vect solve_via_LU(const Vect&);
  Vect solve_via_Cholesky(const Vect&);

  DenseMatrix& operator=(const DenseMatrix&) = default;
  DenseMatrix& operator=(DenseMatrix&&) = default;
  double& operator()(int, int);
  double operator()(int, int) const;



  friend std::ostream& operator<<(std::ostream&, const DenseMatrix&);
  friend DenseMatrix operator*(double, const DenseMatrix&);
  friend DenseMatrix operator+(const DenseMatrix&, const DenseMatrix&);
  friend DenseMatrix operator-(const DenseMatrix&, const DenseMatrix&);
  friend DenseMatrix operator*(const DenseMatrix&, const DenseMatrix&);
  friend Vect operator*(const DenseMatrix&, const Vect&);

};


template<class U, class V, class W>
class tuple3;

template<class U, class V, class W>
bool operator<(const tuple3<U, V, W>&, const tuple3<U, V, W>&);

template<class U, class V, class W>
bool operator>(const tuple3<U, V, W>&, const tuple3<U, V, W>&);

template<class U, class V, class W>
bool operator==(const tuple3<U, V, W>&, const tuple3<U, V, W>&);

template<class U, class V, class W>
class tuple3{
public:
  U first;
  V second;
  W third;

  friend bool operator< <>(const tuple3&, const tuple3&);
  friend bool operator> <>(const tuple3&, const tuple3&);
  friend bool operator== <>(const tuple3&, const tuple3&);
};

class COO : public Matrix{

  std::vector<tuple3<int, int, double> > data;

public:

  COO(const DenseMatrix&);

  double& operator()(int, int) override;
  double operator()(int, int) const override;

  friend Vect operator*(const COO&, const Vect&);
};






#endif
