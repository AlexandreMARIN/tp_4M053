#ifndef MATRIX_HPP
#define MATRIX_HPP

#include <vector>
#include <iostream>


#include "Vect.hpp"


class Matrix;
class COO;
class CSR;


class Matrix{
public:
  Matrix() = default;
  Matrix(const Matrix&) = default;
  Matrix(Matrix&&) = default;
  Matrix(int);

  Matrix& operator=(const Matrix&) = default;
  Matrix& operator=(Matrix&&) = default;

  virtual double& operator()(int, int) = 0;
  virtual double operator()(int, int) const = 0;

  int size() const;

protected:
  int N_ = 0;//number of lines and number of columns

};

Vect operator*(const Matrix&, const Vect&);


class DenseMatrix : public Matrix{
private:
  enum state_ {NORMAL, LU, CHOLESKY} state = NORMAL;


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

template<class U, class V, class W>
bool operator<(const tuple3<U, V, W>& t1, const tuple3<U, V, W>& t2){
  if(t1.first<t2.fist || (t1.first==t2.first && (t1.second<t2.second || (t1.second==t2.second && t1.third<t2.third) ) ) ){
    return true;
  }
  return false;
}

template<class U, class V, class W>
bool operator>(const tuple3<U, V, W>& t1, const tuple3<U, V, W>& t2){
  if(t1.first>t2.fist || (t1.first==t2.first && (t1.second>t2.second || (t1.second==t2.second && t1.third>t2.third) ) ) ){
    return true;
  }
  return false;
}

template<class U, class V, class W>
bool operator==(const tuple3<U, V, W>& t1, const tuple3<U, V, W>& t2){
  if(t1.first==t2.fist && t1.second==t2.second && t1.third==t2.third){
    return true;
  }
  return false;
}


class COO : public Matrix{

  std::vector<tuple3<int, int, double> > data;

public:

  COO(const DenseMatrix&);

  static COO Laplacian(int n);

  double& operator()(int, int) override;
  double operator()(int, int) const override;

  void display() const;
  int get_nnz() const;//number of non-zeros

  friend class CSR;
  friend Vect operator*(const COO&, const Vect&);
};


class CSR : public Matrix{

  std::vector<int> row, col;
  std::vector<double> val;

public:

  CSR(const COO&);


  double& operator()(int, int) override;
  double operator()(int, int) const override;


  friend Vect operator*(const CSR&, const Vect&);
};




#endif
