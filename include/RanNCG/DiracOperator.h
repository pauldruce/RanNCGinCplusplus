//
//  DiracOperator.h
//  RandomNCG
//
//  Created by Paul Druce on 25/12/2021.
//

#ifndef DiracOperator_h
#define DiracOperator_h

#include <map>
#include <armadillo>
#include "Clifford.h"

class DiracOperator
{
private:
   Clifford clifford_module;
   // The size of the matrices inside the commutators and anti-commutators.
   int matrix_size;
   // The size of the dirac operator when fully assembled.
   int dirac_dim;

   // Matrix representation of the Dirac Operator;
   arma::cx_mat dirac_op_mat;

   // Off products of the gamma matrices enter the dirac operator. This function generates them.
   void generate_odd_gamma_products();
   std::vector<arma::cx_mat> odd_gamma_products;

   // Pairs of hermitian gamma products and hermitian matrices;
   std::vector<std::pair<arma::cx_mat, arma::cx_mat>> gamma_H_pairs;
   std::vector<std::pair<arma::cx_mat, arma::cx_mat>> gamma_H_pairs_new;

   // Pairs of anti-hermitian gamma products and anti-hermitian matrices;
   std::vector<std::pair<arma::cx_mat, arma::cx_mat>> gamma_L_pairs;
   std::vector<std::pair<arma::cx_mat, arma::cx_mat>> gamma_L_pairs_new;

public:
   DiracOperator(Clifford &clifford_module, int &N);
   void update();
   // Returns a random Dirac operator with entries uniformaly sampled between [-1,1] + i[-1,1].
   void random_dirac(double step_size);

   DiracOperator &operator+=(DiracOperator &D_right);
   void add(arma::cx_mat &delta_D, std::vector<std::pair<arma::cx_mat, arma::cx_mat>> gamma_HL_pairs);

   inline int get_dirac_dim() { return dirac_dim; };
   inline arma::cx_mat &as_matrix() { return dirac_op_mat; };
   inline int get_matrix_size() { return matrix_size; };
   inline std::vector<arma::cx_mat> &get_odd_gamma_products() { return odd_gamma_products; };
   inline std::vector<std::pair<arma::cx_mat, arma::cx_mat>> get_anti_herm_pairs() { return gamma_L_pairs; };
   inline std::vector<std::pair<arma::cx_mat, arma::cx_mat>> get_herm_pairs() { return gamma_H_pairs; };
   // inline std::vector<std::pair<arma::cx_mat, arma::cx_mat>> get_gamma_HL_pairs(){return gamma_HL_pairs;};
};

#endif /* DiracOperator_h */
