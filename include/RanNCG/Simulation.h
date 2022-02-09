//
//  Simulation.hpp
//  RandomNCG
//
//  Created by Paul Druce on 29/12/2021.
//

#ifndef Simulation_h
#define Simulation_h

#include "DiracOperator.h"
#include <armadillo>

class Simulation
{
private:
   int num_moves;
   int accepted_moves;
   double acceptance_rate;
   double action_val;
   double step_size;

   double Action(DiracOperator &dirac);
   void Metropolis();
   DiracOperator D;
   DiracOperator proposed_D;
   double g2, g4;
   double weight;

public:
   Simulation(DiracOperator dirac_operator);
   void set_params(double &g2, double &g4, double &weight);
   double run_simulation(int chain_length = 1000, double step_size = 0.01);
   inline double get_S() { return action_val; };
   inline arma::cx_mat get_dirac_op() { return D.as_matrix(); };
};

#endif /* Simulation_hpp */
