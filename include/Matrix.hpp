#ifndef MATRIX_HPP
#define MATRIX_HPP

#include <vector>
#include <iostream>


#include "Vect.hpp"



class Matrix{
private:
  int N_;//number of lines and number of columns
  std::vector<double> coef_;//entries
public:
  Matrix();
  Matrix(int);
  Matrix(const Matrix&) = default;
  ~Matrix() = default;

  int size() const;
  void display() const;
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
