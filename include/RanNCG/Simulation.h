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
protected:
   // Default Action calculation. This is slow, you can speed up the simulation by overriding this function.
   virtual double Action(DiracOperator &dirac) const;

public:
   Simulation(const DiracOperator &dirac_operator, SimulationData &simData);

   double run_simulation(int chain_length = 1000, double step_size = 0.01, bool record_action = false);

   [[nodiscard]] inline double get_S() const { return action_val; };
   inline arma::cx_mat get_dirac_op() { return D.as_matrix(); };
   void reset_dirac();

protected:
   SimulationData &sim_data;
private:
   int num_moves;
   int accepted_moves;
   double action_val;
   

   void Metropolis();
   DiracOperator D;
   DiracOperator proposed_D;
};

#endif /* Simulation_h */
