
#include <string>
#include <filesystem>
#include <armadillo>

using namespace arma;
namespace fs = std::filesystem;

#define PROJECT_PATH "/Users/pauldruce/Dev/RanNCGinC++"

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
   std::string action_value_output = "";
   std::string dirac_dataset_name = "";
   std::string hdf_filename = "";
   int dirac_id = 0;
   
public:
   void create_action_output_file_if_needed()
   {
      std::ostringstream filename;
      
      auto t = std::time(nullptr);
      auto tm = *std::localtime(&t);
      
      filename << "action_data_" << this->p << "_" << this->q << "_N_" << this->matrix_size << "_g2_" << -1.0 * g2 << std::put_time(&tm, "_%d-%m-%Y") << ".txt";
      action_value_output = filename.str();
   }
   
   void create_dirac_hdf_name()
   {
      if(hdf_filename == "")
      {
         std::ostringstream filename;
         filename << "dirac_matrices_" << this->p << "_" << this->q << "_N_" << this->matrix_size << ".h5";
         hdf_filename = filename.str();
      }
   }
   
   void create_dirac_dataset_name()
   {
      std::ostringstream filename;
      
      auto t = std::time(nullptr);
      auto tm = *std::localtime(&t);
      
      filename << "/dirac_matrices"
      << "_g2_" <<std::fixed <<std::setprecision(3) << -1.0 * g2 << std::put_time(&tm, "_%Y-%m-%d");
      
      filename << "/dirac_" << dirac_id;
      dirac_id++;
      dirac_dataset_name = filename.str();
   }
   
   void create_output_dir_if_needed()
   {
      if (!fs::exists((std::string)PROJECT_PATH + "/output/"))
      {
         fs::create_directory((std::string)PROJECT_PATH + "/output/");
      }
   }
   
   void print_step_size()
   {
      std::ostringstream filename;
      filename << (std::string)(PROJECT_PATH) + "/output/step_sizes_type_" << this->p << "_" << this->q << ".txt";
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
      file  << std::fixed << std::setprecision(3) << this-> g2 <<","
            << this->matrix_size << ","
            << std::fixed << std::setprecision(17) << this->step_size << std::endl;
      file.close();
   }
   
   void print_action_data()
   {
      // std::ostringstream filename;
      
      create_action_output_file_if_needed();
      create_output_dir_if_needed();
      std::string filepath = (std::string)PROJECT_PATH + "/output/" + action_value_output;
      
      std::ofstream file;
      file.precision(17);
      file.open(filepath, ios::out | ios::app);
      if (!file.is_open())
      {
         std::cout << "ERROR: Could not open file " << filepath << std::endl;
         return;
      }
      file << this->action_value << ", " << this->acceptance_rate << std::endl;
      file.close();
   }
   
   void print_dirac_op_data(const cx_mat dirac_op)
   {
      create_dirac_dataset_name();
      create_output_dir_if_needed();
      create_dirac_hdf_name();
      std::string hdf_filepath = (std::string)PROJECT_PATH + "/output/" + hdf_filename;
      dirac_op.save(hdf5_name(hdf_filepath, dirac_dataset_name, hdf5_opts::append));
   }
   
   
};
