#include <iostream>
#include <armadillo>
#include <filesystem>
#include <string>

//#include "Clifford.h"
//#include "DiracOperator.h"
//#include "Simulation.h"

#define PROJECT_PATH "/Users/pauldruce/Dev/RanNCGinC++"

std::vector<arma::cx_mat> LoadDiracOperators(const std::string& filename, const std::string& data_set_name)
{
    std::vector<arma::cx_mat> dirac_opts;
    arma::cx_mat op;
    for(int i = 1; i< 21; i++)
    {
        std::string data_set_name_final = data_set_name + "/dirac_" + std::to_string(i);
        if(op.load(arma::hdf5_name(filename, data_set_name_final))) {
            dirac_opts.push_back(op);
        }
        else
        {
            throw std::runtime_error("Error reading matrix from hdf5 file.");
        }
    }

    return dirac_opts;
}

int main()
{
    // Load dirac operators
    double g2 = 0.0;
    std::string filename = (std::string)(PROJECT_PATH) + "/output/dirac_matrices_1_1_N_5.h5";
    std::stringstream  g2_string;
    g2_string << std::fixed << std::setprecision(3) <<g2;
    std::string data_set_name = "dirac_matrices_g2_" + g2_string.str() + "_2022-02-06";

    std::vector<arma::cx_mat> diracs = LoadDiracOperators(filename, data_set_name);


    // calculate eigenvalues
    std::vector<arma::vec> all_eigs;
    for(int i = 0; i<20; i++) {
        arma::vec eigs = eig_sym(diracs[0]);
        all_eigs.push_back(eigs);
    }

    // take average of eigenvalues with errors. (errors will likely be large for now. )
    std::vector<double> avg_eigs;
    for(int i=0; i<all_eigs[0].size(); i++)
    {
        double avg=0;

        for(auto eigs :all_eigs)
        {
            avg+= eigs[i];
        }

        avg /= all_eigs.size();
        avg_eigs.push_back(avg);
    }


    for(auto e : avg_eigs) {
        std::cout << e << std::endl;
    }

    // export eigenvalues to file.
    // #TODO: Then create python notebook to analyse eigenvalues distribution.
    return -1;
}