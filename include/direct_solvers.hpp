#ifndef DIRECT_SOLVERS_HPP
#define DIRECT_SOLVERS_HPP


#include "Matrix.hpp"

Vect solve_triang_inf(const Matrix&, const Vect&);
Vect solve_triang_sup(const Matrix&, const Vect&);
Vect solve_triang_inf_id(const Matrix&, const Vect&);
Vect solve_triang_sup_id(const Matrix&, const Vect&);




#endif
