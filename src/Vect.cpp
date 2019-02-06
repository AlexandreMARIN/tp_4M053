#include <cmath>
#include "../include/Vect.hpp"

using namespace std;

Vect::Vect() : N_(3), coef_(3){}

Vect::Vect(int N) : N_(N), coef_(N){}

Vect::Vect(initializer_list<double> il) : N_(il.size()), coef_(il.size()){
  for(int i=0;i<N_;i++){
    coef_[i] = il.begin()[i];
  }
}

int Vect::size() const{
  return N_;
}

void Vect::display() const{
  cout << scientific << setprecision(2) << "N_ : " << N_ << "\ncoef_ :\n( ";
  if(coef_.size()){
    cout << setw(10) << coef_[0];
  }
  for(unsigned int i=1;i<coef_.size();i++){
    cout << ", " << setw(10) << coef_[i];
  }
  cout << " )\n\n";
}

void Vect::resize(int NN_){
  N_ = NN_;
  coef_.resize(NN_);
}

void Vect::set_to_zero(){
  for(unsigned int i=0;i<coef_.size();i++){
    coef_[i] = 0.0;
  }
}

void Vect::fill(double alpha){
  for(int i=0;i<N_;i++){
    coef_[i] = alpha;
  }
}

double Vect::norm_infty() const{
  double res = abs(coef_[0]), aux;

  for(int i=1;i<N_;i++){
    aux = abs(coef_[i]);
    if(aux>res){
      res = aux;
    }
  }

  return res;
}

double& Vect::operator()(int i){
  return coef_[i];
}

double Vect::operator()(int i) const{
  return coef_[i];
}

const Vect& Vect::operator+=(const Vect& v){
  for(int i=0;i<N_;i++){
    coef_[i] += v.coef_[i];
  }
  return *this;
}

const Vect& Vect::operator*=(double alpha){
  for(int i=0;i<N_;i++){
    coef_[i] *= alpha;
  }
  return *this;
}


ostream& operator<<(ostream& os, const Vect& v){
  os << "N_ : " << v.N_ << "\ncoef_ :\n( ";
  if(v.coef_.size()){
    os << v.coef_[0];
  }
  for(unsigned int i=1;i<v.coef_.size();i++){
    os << ", " << v.coef_[i];
  }
  os << " )\n\n";

  return os;
}

Vect operator+(const Vect& v1, const Vect& v2){

  Vect sum(v1.N_);
  for(int i=0;i<v1.N_;i++){
    sum.coef_[i] = v1.coef_[i] + v2.coef_[i];
  }
  return sum;

}

Vect operator-(const Vect& v1, const Vect& v2){

  Vect sum(v1.N_);
  for(int i=0;i<v1.N_;i++){
    sum.coef_[i] = v1.coef_[i] - v2.coef_[i];
  }
  return sum;

}

Vect operator*(double alpha, const Vect& v){
  Vect res(v.N_);
  for(unsigned int i=0;i<v.coef_.size();i++){
    res.coef_[i] = alpha*v.coef_[i];
  }
  return res;
}
