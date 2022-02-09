#pragma once
#include "Clifford.h"
#include <armadillo>

class Pauli : protected Clifford
{
   arma::cx_mat x, y, z;
   Pauli();
   
public:
   void introduce() override;
};
