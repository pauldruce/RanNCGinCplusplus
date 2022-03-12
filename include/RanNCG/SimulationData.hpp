//
//  SimulationData.hpp
//  RandomNCG
//
//  Created by Paul Druce on 29/12/2021.
//

#ifndef SimulationData_hpp
#define SimulationData_hpp

#include <string>
#include <filesystem>
#include <armadillo>

using namespace arma;
namespace fs = std::filesystem;

struct SimulationData
{
   int p;
   int q;
   double g2;
   double g4;
   int matrix_size;
   double step_size = 0.0;
   double action_value = 0.0;
   double acceptance_rate = 0.0;

   std::string project_path;
   std::string output_path;

   std::string action_folder;
   std::string action_values_filename;

   int dirac_id = 0;
   std::string dirac_folder;
   std::string dirac_dataset_name;
   std::string dirac_hdf_filename;

   int eigenvalues_id = 0;
   std::string eigenvalues_folder;
   std::string eigenvalues_dataset_name;
   std::string eigenvalues_hdf_filename;

public:
   void create_output_dir()
   {
      if(!fs::exists(project_path+"/output/"))
      {
         fs::create_directory(project_path+"/output/");
      }
      std::ostringstream ss;
      ss << "/output/"
         << "Simulation_" << this->p << "_" << this->q << "_N_" << this->matrix_size << "/";
      output_path = project_path + ss.str();
      if (!fs::exists(output_path))
      {
         fs::create_directory(output_path);
      }
   }

   void create_action_folder()
   {
      if (output_path.empty())
      {
         create_output_dir();
      }

      std::ostringstream ss;
      ss << "action_data" << this->p << "_" << this->q << "_N_" << this->matrix_size << "/";
      action_folder = output_path + ss.str();
      create_output_dir();
      if (!fs::exists(action_folder))
      {
         fs::create_directory(action_folder);
      }
   }
   void create_dirac_dir()
   {
      if (output_path.empty())
      {
         create_output_dir();
      }

      std::ostringstream ss;
      ss << "dirac_matrices_" << this->p << "_" << this->q << "_N_" << this->matrix_size << "/";
      dirac_folder = output_path + ss.str();
      if (!fs::exists(dirac_folder))
      {
         fs::create_directory(dirac_folder);
      }
   }

   void create_eigenvalues_dir()
   {
      if (output_path.empty())
      {
         create_output_dir();
      }

      std::ostringstream ss;
      ss << "eigenvalues_" << this->p << "_" << this->q << "_N_" << this->matrix_size << "/";
      eigenvalues_folder = output_path + ss.str();
      if (!fs::exists(eigenvalues_folder))
      {
         fs::create_directory(eigenvalues_folder);
      }
   }

   void create_action_filename()
   {
      std::ostringstream ss;
      auto t = std::time(nullptr);
      auto tm = *std::localtime(&t);

      ss << std::fixed << std::setprecision(3);
      ss << "action_data_"
         << this->p << "_" << this->q << "_"
         << "N_" << this->matrix_size
         << "_g2_" << -1.0 * g2
         << std::put_time(&tm, "_%Y-%m-%d")
         << ".csv";

      action_values_filename = ss.str();
   }

   std::string create_dirac_hdf_filename()
   {
      std::ostringstream ss;
      ss << "Dirac_Matrices_"
         << this->p << "_" << this->q
         << "_N_" << this->matrix_size
         << std::fixed << std::setprecision(3) << "_g2_" << -1.0 * g2
         << ".h5";
      dirac_hdf_filename = ss.str();
      return dirac_hdf_filename;
   }

   std::string create_eigenvalues_hdf_filename()
   {
      std::ostringstream ss;
      ss << "Eigenvalues_"
         << this->p << "_" << this->q
         << "_N_" << this->matrix_size
         << std::fixed << std::setprecision(3) << "_g2_" << -1.0 * g2
         << ".h5";
      eigenvalues_hdf_filename = ss.str();
      return eigenvalues_hdf_filename;
   }

   void create_dirac_dataset()
   {
      std::ostringstream ss;
      ss << "/dirac_" << dirac_id;
      dirac_id++;
      dirac_dataset_name = ss.str();
   }

   void create_eigenvalues_dataset()
   {
      std::ostringstream ss;
      ss << "/eigenvalues_" << eigenvalues_id;
      eigenvalues_id++;
      eigenvalues_dataset_name = ss.str();
   }

   void print_step_size()
   {
      create_output_dir();
      std::ostringstream filename;
      filename << output_path + "/step_sizes_type_" << this->p << "_" << this->q << ".txt";
      // a+ means pointer begins at beginning of file BUT when written to, it will always be at the end of the file.
      // Using a+ because I may want to search and update values.
      std::fstream file;
      file.precision(17);
      file.open(filename.str(), ios::in | ios::out | ios::app);
      if (!file.is_open())
      {
         std::cout << "ERROR: Failed to open file with name '" << filename.str() << "'\n";
         return;
      }
      file << std::fixed << std::setprecision(3) << this->g2 << ","
           << this->matrix_size << ","
           << std::fixed << std::setprecision(17) << this->step_size << std::endl;
      file.close();
   }

   void print_action_data()
   {
      create_output_dir();
      create_action_folder();
      create_action_filename();
      std::string filepath = action_folder + action_values_filename;
      std::fstream file;
      file.open(filepath, ios::in | ios::app);
      if (!file.is_open())
      {
         std::cout << "ERROR: Failed to open file with name '" << filepath << "'\n";
         return;
      }
      file << std::fixed << std::setprecision(17)
           << this->action_value << ", " << acceptance_rate << std::endl;
      file.close();
   }

   void print_dirac_op_data(const cx_mat &dirac_op)
   {
      create_output_dir();
      create_dirac_dir();
      create_dirac_hdf_filename();
      create_dirac_dataset();

      std::string hdf_filepath = dirac_folder + dirac_hdf_filename;
      dirac_op.save(hdf5_name(hdf_filepath, dirac_dataset_name, hdf5_opts::append));
   }

   void save_eigenvalues(const cx_mat &dirac_op)
   {
      create_output_dir();
      create_eigenvalues_dir();
      create_eigenvalues_hdf_filename();
      create_eigenvalues_dataset();

      std::string hdf_filepath = eigenvalues_folder + eigenvalues_hdf_filename;
      // get eigenvalue
      vec eigs;
      if(!dirac_op.is_hermitian(1e-10))
      {
         auto hermit_dirac = dirac_op;
         hermit_dirac.clean(1e-10);
         eigs = eig_sym(hermit_dirac);
      }
      else
      {
         eigs = eig_sym(dirac_op);
      }
      // save eigs
      eigs.save(hdf5_name(hdf_filepath, eigenvalues_dataset_name, hdf5_opts::append));
   }
};

#endif /* SimulationData_hpp */
