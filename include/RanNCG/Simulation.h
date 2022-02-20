//
//  Simulation.hpp
//  RandomNCG
//
//  Created by Paul Druce on 29/12/2021.
//

#ifndef Simulation_h
#define Simulation_h

#include "DiracOperator.h"
#include "SimulationData.hpp"
#include <armadillo>

class Simulation
{
private:
   int num_moves{};
   int accepted_moves{};
   double acceptance_rate{};
   double action_val{};

   SimulationData& sim_data;
//   double step_size;
//   double g2, g4;
//   std::string hdf_file;
//   std::string action_dataset;

   double Action(DiracOperator &dirac);
   void Metropolis();
   DiracOperator D;
   DiracOperator proposed_D;

public:
   Simulation(DiracOperator dirac_operator, SimulationData &simData);
   double run_simulation(int chain_length = 1000, double step_size = 0.01, bool record_action = false);
   inline double get_S() { return action_val; };
   inline arma::cx_mat get_dirac_op() { return D.as_matrix(); };

};

#endif /* Simulation_hpp */
