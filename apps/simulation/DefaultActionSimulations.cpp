#include <iostream>
#include <armadillo>
#include <string>
#include <cassert>

#include "Clifford.h"
#include "DiracOperator.h"
#include "Simulation.h"
#include "SimulationData.hpp"

using namespace arma;
namespace fs = std::filesystem;

int main()
{
   // Set the project path
   std::string project_path = "/Users/pauldruce/Dev/RanNCGinC++";
   // Initialize the random number generator so that it is actually random.
   arma_rng::set_seed_random();

   // Set Clifford type parameters for simulation;
   int p = 1;
   int q = 3;
   Clifford cliff(p, q);

   // Set up the action parameters
   double g4 = 1.0;
   double g2;
   double g2_start = -0.00;
   double g2_end = -4.0;
   double g2_step = -0.04;

   // Set simulation parameters for initial stage of finding an appropriate step size.
   int chain_length = 200;
   int num_runs_before_reset = 1000;
   bool record_action = false;
   double step_size = 0.0243381;

   // Set size of H and L matrices to use.
   int matrix_size;
   for (matrix_size = 4; matrix_size < 11; matrix_size++)
   {

      for (g2 = g2_start; g2 > g2_end; g2 += g2_step)
      {
         SimulationData sim_data;
         sim_data.project_path = project_path;
         sim_data.matrix_size = matrix_size;
         sim_data.p = p;
         sim_data.q = q;
         sim_data.g2 = g2;
         sim_data.g4 = g4;
         sim_data.step_size = step_size;

         DiracOperator dirac_operator(cliff, sim_data.matrix_size);
         Simulation simulation(dirac_operator, sim_data);

         cout << "Simulation parameters: (g2, g4, N) = (" << g2 << "," << g4 << "," << sim_data.matrix_size << ")"
              << std::endl;

         int num_times_acceptance_ratio_is_okay = 0;

         double acceptance_rate_tol = 0.01;
         double acceptance_rate_upper_bound = 0.5 * (1.0 + acceptance_rate_tol);
         double acceptance_rate_lower_bound = 0.5 * (1.0 - acceptance_rate_tol);
         int num_tempering_runs = 0;

         while (sim_data.acceptance_rate > acceptance_rate_upper_bound ||
                sim_data.acceptance_rate < acceptance_rate_lower_bound ||
                num_times_acceptance_ratio_is_okay < 10)
         {
            if (num_tempering_runs % num_runs_before_reset == 0 or sim_data.step_size < 1e-12)
            {
               // This occurs if there have been a large number of runs (~200,000) runs or if the step size gets too small.
               sim_data.step_size = step_size;
               simulation.reset_dirac();
            }
            // Make sure step_size is positive
            sim_data.step_size = std::abs(sim_data.step_size);

            assert(sim_data.step_size > 1e-12);

            sim_data.acceptance_rate = simulation.run_simulation(chain_length, sim_data.step_size, record_action);
            sim_data.action_value = simulation.get_S();
            if (num_tempering_runs % 100 == 0)
            {
               std::cout << "tempering_run = " << num_tempering_runs
                         << ", step_size = " << sim_data.step_size
                         << ", action_val = " << sim_data.action_value
                         << ", acceptance ratio = " << sim_data.acceptance_rate << std::endl;
            }

            if (sim_data.acceptance_rate > acceptance_rate_upper_bound or sim_data.acceptance_rate < acceptance_rate_lower_bound)
            {
               sim_data.step_size += ((sim_data.acceptance_rate - 0.5) * sim_data.step_size);
            }
            else
            {
               num_times_acceptance_ratio_is_okay++;
            }

            num_tempering_runs++;
         }

         std::cout << "System has found appropriate step size for configuration: "
                   << sim_data.step_size << std::endl;

         sim_data.print_step_size();
         std::cout << "Running the simulation for 5000 moves to burn in." << std::endl;
         sim_data.acceptance_rate = simulation.run_simulation(5000, sim_data.step_size, true);

         int num_runs = 2000;
         int num_moves_per_run = 4 * sim_data.matrix_size;
         std::cout << "Running the simulation for " << num_runs << " many runs, each with " << num_moves_per_run
                   << " moves." << std::endl;
         for (int i = 0; i < num_runs; i++)
         {
            sim_data.acceptance_rate = simulation.run_simulation(num_moves_per_run, sim_data.step_size, true);
            // sim_data.print_dirac_op_data(simulation.get_dirac_op());
            sim_data.save_eigenvalues(simulation.get_dirac_op());
         }

         step_size = sim_data.step_size;
      }
   }

   return 0;
}
