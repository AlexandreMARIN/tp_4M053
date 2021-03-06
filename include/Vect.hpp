#ifndef VECT_HPP
#define VECT_HPP

#include <vector>
#include <iostream>
#include <stdexcept>
#include <iomanip>


class Vect{
private:
  int N_;
  std::vector<double> coef_;
public:
  Vect();
  Vect(int);
  Vect(const Vect&) = default;
  Vect(Vect&&) = default;
  Vect(std::initializer_list<double>);
  ~Vect() = default;

  int size() const;
  void display() const;

  void resize(int);
  void set_to_zero();
  void fill(double);

  double norm_infty() const;

  Vect& operator=(const Vect&) = default;
  Vect& operator=(Vect&&) = default;
  double& operator()(int);
  double operator()(int) const;
  const Vect& operator+=(const Vect&);
  const Vect& operator-=(const Vect&);
  const Vect& operator*=(double);


  friend std::ostream& operator<<(std::ostream&, const Vect&);
  friend Vect operator+(const Vect&, const Vect&);
  friend Vect operator-(const Vect&, const Vect&);
  friend Vect operator*(double, const Vect&);
  friend double operator,(const Vect&, const Vect&);
};




#endif
