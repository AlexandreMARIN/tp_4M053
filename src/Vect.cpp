#include "../include/Vect.hpp"

using namespace std;

Vect::Vect() : N_(3), coef_(3){}

Vect::Vect(int N) : N_(N), coef_(N){}

/*Vect::Vect(const Vect& v) : N_(v.N_), coef_(v){

  }*/

int Vect::size() const{
  return N_;
}

void Vect::display() const{
  cout << "N_ : " << N_ << "\ncoef_ :\n( ";
  if(coef_.size()){
    cout << coef_[0];
  }
  for(unsigned int i=1;i<coef_.size();i++){
    cout << ", " << coef_[i];
  }
  cout << " )\n\n";
}

double& Vect::operator()(int i){
  return coef_[i];
}

double Vect::operator()(int i) const{
  return coef_[i];
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

Vect operator*(double alpha, const Vect& v){
  Vect res(v.N_);
  for(unsigned int i=0;i<v.coef_.size();i++){
    res.coef_[i] = alpha*v.coef_[i];
  }
  return res;
}
