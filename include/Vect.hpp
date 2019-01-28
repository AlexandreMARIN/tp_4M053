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

  Vect& operator=(const Vect&) = default;
  Vect& operator=(Vect&&) = default;
  double& operator()(int);
  double operator()(int) const;



  friend std::ostream& operator<<(std::ostream&, const Vect&);
  friend Vect operator*(double, const Vect&);

};












#endif
