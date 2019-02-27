#ifndef MATRIX_HPP
#define MATRIX_HPP

#include <vector>
#include <iostream>


#include "Vect.hpp"



class Matrix{
private:
  int N_;//number of lines and number of columns
  std::vector<double> coef_;//entries
  bool LU_ = false;
  std::vector<int> perm_;
public:
  Matrix();
  Matrix(int);
  Matrix(const Matrix&) = default;
  Matrix(Matrix&&) = default;
  ~Matrix() = default;

  int size() const;
  void display() const;
  void display_LU() const;

  void diag(double, int=0);
  void decomp_LU();
  Vect solve_via_LU(const Vect&);

  Matrix& operator=(const Matrix&) = default;
  Matrix& operator=(Matrix&&) = default;
  double& operator()(int, int);
  double operator()(int, int) const;



  friend std::ostream& operator<<(std::ostream&, const Matrix&);
  friend Matrix operator*(double, const Matrix&);
  friend Matrix operator+(const Matrix&, const Matrix&);
  friend Matrix operator-(const Matrix&, const Matrix&);
  friend Matrix operator*(const Matrix&, const Matrix&);
  friend Vect operator*(const Matrix&, const Vect&);

};












#endif
