//
//  Simulation.cpp
//  RandomNCG
//
//  Created by Paul Druce on 29/12/2021.
//

#include "Simulation.h"
#include "SimulationData.hpp"

using namespace arma;
using namespace std;

Simulation::Simulation(const DiracOperator & dirac_operator, SimulationData &simData)
      : D(dirac_operator),
        proposed_D(dirac_operator), sim_data(simData)
{
//   this->g2 = sim_data.g2;
//   this->g4 = sim_data.g4;
//   this->step_size = sim_data.step_size;
}

double Simulation::run_simulation(int chain_length, double step_size, bool record_action)
{
   this->sim_data.step_size = step_size;

   // Reset the acceptance rate values;
   this->accepted_moves = 0;
   this->num_moves = 0;
//   this->acceptance_rate = 0;

   // Get the initial action value for comparison later.
   action_val = Action(this->D);

   for (int i = 0; i < chain_length; i++)
   {
      Metropolis();
      sim_data.acceptance_rate = (double) this->accepted_moves / (double) this->num_moves;
      if (record_action)
      {
         sim_data.action_value = this->action_val;
         sim_data.print_action_data();
      }
   }
   // std::cout << "Completed Metropolis" << std::endl;
   //   this->acceptance_rate = (double) this->accepted_moves / (double) this->num_moves;
   // std::cout << acceptance_rate << std::endl;
   //   return this->acceptance_rate;
   return this->sim_data.acceptance_rate;
}

void Simulation::Metropolis()
{
   // Updates the Dirac operator according to the Metropolis-Hastings algorithm
   // This function inputs the old Dirac operator, the action coupling constants, weightA
   // which dictates the steps size between the old and new Dirac operators, and a variable to
   // calculate the acceptance_rate. This function could be done smoother I'm sure. But it seems to work.

   // Propose a new Dirac operator;
   proposed_D.random_dirac(sim_data.step_size);
   double S_new = Action(proposed_D);
   // I think of rand_press as the random pressure which pushes you out of sufficiently small local minima.
   double delta_S = this->action_val - S_new;
   double rand_press = exp(delta_S);
   // p is the random value between 0 and 1 used in Metropolis-Hastings algorithm
   double p = randu();

   // This is my understanding of the Metropolis-Hastings algorithm
   if (S_new < this->action_val or rand_press > p)
   {
      D = proposed_D;
      this->action_val = S_new;
      this->accepted_moves += 1;
      this->num_moves += 1;
   } else
   {
      // Don't update the dirac operator.
      proposed_D = D;
      this->num_moves += 1;
   }
}

double Simulation::Action(DiracOperator &dirac) const
{
   // Currently, using D matrix. Need to migrate to H and L matrices.
   cx_mat D_mat = dirac.as_matrix();
   cx_mat D2 = D_mat * D_mat;

   cx_mat D4 = D2*D2;

   cx_double action = sim_data.g2 * trace(D2) + sim_data.g4 * trace(D4);
   return action.real();
}

void Simulation::reset_dirac()
{
   this->D.reset_dirac();
   this->action_val = Action(D);
}
