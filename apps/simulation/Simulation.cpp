#include <iostream>
#include <armadillo>
#include <filesystem>
#include <string>
#include <fstream>

#include "Clifford.h"
#include "DiracOperator.h"
#include "Simulation.h"

#include "SimulationData.hpp"

// #TODO: Need to make the outfile have information about the g2 values.

// using namespace std;
using namespace arma;
namespace fs = std::filesystem;

int main()
{
   // Initialize the random number generator so that it is actually random.
   arma_rng::set_seed_random();

   double g4 = 1.0;
   double g2 = 0.0;
   double g2_start = -0.0;
   double g2_end = -4.0;
   double num_steps = 100;
//   double g2_step = (g2_end - g2_start) / num_steps;
   double g2_step = -0.04;
   int chain_length = 500;
   int matrix_size= 9;
//   for (int matrix_size = 8; matrix_size < 15; matrix_size++)
//   {
      for (g2 = g2_start; g2 > g2_end; g2 += g2_step)
      {

         SimulationData sim_data;
         sim_data.matrix_size = matrix_size;
         sim_data.p = 1;
         sim_data.q = 1;
         sim_data.g2 = g2;
         sim_data.g4 = g4;
//         sim_data.step_size = 1.0/5.0/std::pow(sim_data.matrix_size, 1.5);
         sim_data.step_size = 0.005;
         sim_data.create_output_dir_if_needed();
         sim_data.create_dirac_hdf_name();
         sim_data.create_dirac_dataset_name();
         sim_data.create_action_output_file_if_needed();
         
         Clifford cliff(sim_data.p, sim_data.q);
         cliff.introduce();
         DiracOperator dirac_operator(cliff, sim_data.matrix_size);
         Simulation simulation(dirac_operator);

         // Run the simulation for a few runs, check the step size.

         cout << "Simulation parameters: (g2, g4, N) = (" << g2 << "," << g4 << "," << sim_data.matrix_size << ")" << std::endl;
         simulation.set_params(g2, g4, sim_data.step_size);
         double step_size_diff;
         int num_times_acceptance_ratio_is_okay = 0;

         double acceptance_rate_tol = 0.01;
         double acceptance_rate_upper_bound = 0.5 * (1.0 + acceptance_rate_tol);
         double acceptance_rate_lower_bound = 0.5 * (1.0 - acceptance_rate_tol);
         while (sim_data.acceptance_rate > acceptance_rate_upper_bound || sim_data.acceptance_rate < acceptance_rate_lower_bound || num_times_acceptance_ratio_is_okay < 10)
         {
            // std::cout << sim_data.step_size << " " << sim_data.acceptance_rate << std::endl;
            step_size_diff = 0.1 * sim_data.step_size + 0.0001;
            
            // Make sure step_size is positive
            sim_data.step_size = std::abs(sim_data.step_size);
            
            sim_data.acceptance_rate = simulation.run_simulation(chain_length, sim_data.step_size);
            sim_data.action_value = simulation.get_S();
            if (sim_data.acceptance_rate > acceptance_rate_upper_bound)
            {
               sim_data.step_size += step_size_diff;
            }
            else if (sim_data.acceptance_rate < acceptance_rate_lower_bound)
            {
               sim_data.step_size -= step_size_diff;
            }
            else
            {
               num_times_acceptance_ratio_is_okay++;
               // std::cout << sim_data.step_size << " " << sim_data.acceptance_rate << std::endl;
            }
         }
         printf("System has found appropriate step size for configuration.\n");
         printf("Outputting step size value for configuration.\n");
         sim_data.print_step_size();
         printf("Running the simulation for 5000 moves to burn in\n");
         sim_data.acceptance_rate = simulation.run_simulation(5000, sim_data.step_size);

         int num_runs = 20;
         int num_moves_per_run = 4 * sim_data.matrix_size;
         std::cout << "Running the simulation for " << num_runs << " many runs, each with " << num_runs*num_moves_per_run << " moves." << std::endl;
         for (int i = 0; i < num_runs; i++)
         {
            sim_data.acceptance_rate = simulation.run_simulation(num_runs * num_moves_per_run, sim_data.step_size);
            sim_data.print_action_data();
            sim_data.print_dirac_op_data(simulation.get_dirac_op());
         }
      }
//   }
   return 0;
}
