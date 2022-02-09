#pragma once

#include "Clifford.h"
#include <armadillo>

class DiracMatrices : protected Clifford
{
private:
   arma::cx_mat gamma0, gamma1, gamma2, gamma3;
   DiracMatrices();

   std::vector<arma::cx_mat> dirac_gens;
   arma::cx_mat dirac_chirality;

   std::vector<arma::cx_mat> chiral_gens;
   arma::cx_mat chiral_chirality;

   std::vector<arma::cx_mat> majorana_gens;
   arma::cx_mat majorana_chirality;

public:
   std::vector<arma::cx_mat> get_dirac_gens() { return dirac_gens; };
   arma::cx_mat get_dirac_chirality() { return dirac_chirality; };

   std::vector<arma::cx_mat> get_chiral_gens() { return chiral_gens; };
   std::vector<arma::cx_mat> get_majorana_gens() { return majorana_gens; };

   void introduce() override;
};
