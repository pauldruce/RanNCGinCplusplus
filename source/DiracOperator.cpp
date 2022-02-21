//
//  DiracOperator.cpp
//  RandomNCG
//
//  Created by Paul Druce on 25/12/2021.
//

#include "DiracOperator.h"
#include "RandomMatrixUtils.h"
#include "cppitertools/combinations.hpp"

using namespace arma;
using namespace std;

DiracOperator::DiracOperator(Clifford &clifford_module, int &N)
    : clifford_module(clifford_module), matrix_size(N)
{
   this->dirac_dim = N * N * clifford_module.get_matrix_size();

   generate_odd_gamma_products();

   // Initialise the gamma-H pairs and gamma-L pairs.
   for (auto prod : odd_gamma_products)
   {
      cx_mat zero_mat(matrix_size, matrix_size, fill::zeros);
      if (prod.is_hermitian()) // Hermitian
      {
         gamma_H_pairs.emplace_back(prod, zero_mat);
      }
      else if ((cx_double(0, 1) * prod).is_hermitian()) // Anti-hermitian
      {
         gamma_L_pairs.emplace_back(prod, zero_mat);
      }
      else
      {
         throw std::logic_error("ERROR: Odd gamma product is neither hermitian or anti-hermitian.");
      }
   }
   // Pairs to hold the proposed new gamma-H and gamma-L pairs when making Monte Carlo move.
   gamma_H_pairs_new = gamma_H_pairs;
   gamma_L_pairs_new = gamma_L_pairs;

   // Set the initial matrix representation of the Dirac to be zero.
   this->dirac_op_mat = zeros<cx_mat>(dirac_dim, dirac_dim);
   random_dirac(1.0);

//this->dirac_op_mat = eye<cx_mat>(dirac_dim, dirac_dim);
}

void DiracOperator::random_dirac(double step_size)
{
   dirac_op_mat.zeros();
   auto herm_pairs = get_herm_pairs();
   auto anti_herm_pairs = get_anti_herm_pairs();

   // Make monte carlo move for the Hermitian matrices.
   for (auto &pair : herm_pairs)
   {

      cx_mat temp = step_size * random_complex(matrix_size);
      cx_mat H_i = temp + temp.t();
      assert(H_i.is_hermitian());
      pair.second += H_i;

      cx_mat anti_comm = anticomm(H_i);
      cx_mat D_temp = kron(pair.first, anti_comm);
      dirac_op_mat += D_temp;
   }

   // Make monte carlo move for the anti-hermitian matrices.
   for (auto &pair : anti_herm_pairs)
   {
      cx_mat temp = step_size * random_complex(matrix_size);
      cx_mat L_i = temp - temp.t();
      assert((cx_double(0, 1) * L_i).is_hermitian());
      // #todo: check if I need to somehow store the new L_i and H_i matrices for a while before adding them.
      pair.second += L_i;
      cx_mat commutator = comm(L_i);
      cx_mat D_temp = kron(pair.first, commutator);
      dirac_op_mat += D_temp;
   }
   assert(dirac_op_mat.is_hermitian(1e-9));
}

DiracOperator &DiracOperator::operator+=(DiracOperator &D_right)
{

   if (this->matrix_size != D_right.matrix_size)
   {
      throw "ERROR: Matrix sizes do not match.";
   }

   if (this->clifford_module.type != D_right.clifford_module.type)
   {
      throw "ERROR: Dirac Operators of different types!";
   }

   auto lhs_herm_pairs = this->get_herm_pairs();
   auto rhs_herm_pairs = D_right.get_herm_pairs();

   for (int i = 0; i < lhs_herm_pairs.size(); i++)
   {
      lhs_herm_pairs[i].second += rhs_herm_pairs[i].second;
   }

   auto lhs_anti_herm_pairs = this->get_anti_herm_pairs();
   auto rhs_anti_herm_pairs = D_right.get_anti_herm_pairs();
   for (int i = 0; i < lhs_herm_pairs.size(); i++)
   {
      lhs_anti_herm_pairs[i].second += rhs_anti_herm_pairs[i].second;
   }

   this->dirac_op_mat += D_right.dirac_op_mat;
   return *this;
}

void DiracOperator::generate_odd_gamma_products()
{
   /* Produce all the odd gamma products

    The function inputs the cliff object (which contains the
    generators) for a given Clifford module, it then extracts the
    uses these to calculate the necessary products
    */
   // Number of odd products (we are doing n choose r for all odd
   // r) number = int(2**k)
   // The plus one is to ensure that if there are an odd num of
   // generators, then the prod of them all is included for i in

   for (int i = 1; i < clifford_module.get_number_of_generators() + 1; i += 2)
   {
      for (auto &&comb : iter::combinations(clifford_module.generators, i))
      {
         if (comb.size() > 1)
         {
            cx_mat product = eye<cx_mat>(clifford_module.get_matrix_size(), clifford_module.get_matrix_size());
            for (auto x : comb)
            {
               product *= x;
            }
            odd_gamma_products.push_back(product);
         }
         else
         {
            odd_gamma_products.push_back(comb[0]);
         }
      }
   }
}
