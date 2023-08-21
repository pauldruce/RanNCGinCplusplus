#ifndef TYPE11SIMULATION_H
#define TYPE11SIMULATION_H

#include "Simulation.h"
#include <armadillo>

class Type11Simulation : public Simulation
{
public:
    Type11Simulation(const DiracOperator &dirac_operator, SimulationData &simData)
        : Simulation(dirac_operator, simData)
    {
        std::pair<int, int> type = dirac_operator.GetType();
        if (type.first != 1 || type.second != 1)
        {
            throw std::runtime_error("Dirac operator is not of the correct type.");
        }
    }

    // Instead of calculating Tr(D^4 + g2 D^2), directly, which involves an (2*N^2) x (2*N^2) matrix for D.
    // We can calculate this action in terms of the subcomponents of D, H and L, which are N x N matrices.
    // This reduces the complexity significantly for larger matrices.
    // See papers by J Barrett and L Glaser on "Monte Carlo simulations of random non-commutative geometries" for details.
    double Action(DiracOperator &dirac_op) const override
    {
        arma::cx_double trH, trH2, trH3, trH4;
        arma::cx_double trL, trL2, trL3, trL4;
        arma::cx_double trHL, trHLHL, trH2L2, trH2L, trL2H;

        auto &all_herm = dirac_op.get_herm_pairs();
        if (all_herm.size() != 1)
            throw std::runtime_error("Incorrect number of Hermitian matrices for type 11 Dirac operator");

        auto &all_anti_herm = dirac_op.get_anti_herm_pairs();
        if (all_anti_herm.size() != 1)
            throw std::runtime_error("Incorrect number of anti-Hermitian matrices for type 11 Dirac operator");

        arma::cx_mat H = all_herm[0].second;
        arma::cx_mat L = all_anti_herm[0].second;

        // Create all the products necessary.
        arma::cx_mat H2 = H * H, H3 = H * H2, H4 = H * H3;
        arma::cx_mat L2 = L * L, L3 = L * L2, L4 = L * L3;
        arma::cx_mat HL = H * L;
        arma::cx_mat HLHL = HL * HL;
        arma::cx_mat H2L2 = H2 * L2;
        arma::cx_mat H2L = H2 * L;
        arma::cx_mat L2H = L2 * H;

        // Calculate the traces necessary.
        arma::cx_double trD2, trD4;

        trH = trace(H);
        trH2 = trace(H2);
        trH3 = trace(H3);
        trH4 = trace(H4);

        trL = trace(L);
        trL2 = trace(L2);
        trL3 = trace(L3);
        trL4 = trace(L4);

        trHL = trace(HL);
        trHLHL = trace(HLHL);
        trH2L2 = trace(H2L2);
        trH2L = trace(H2L);
        trL2H = trace(L2H);

        trD2 = 4.0 * (double)
              H.size() * (trH2 - trL2) + 4.0 * (trH * trH + trL * trL);

        trD4 = 4.0 * (double)H.size() * (trH4 + trL4 - (4.0 * trH2L2) + (2.0 * trHLHL)) + 16.0 * (trH * (trH3 - trL2H) + trL * (-trL3 + trH2L) + trHL * trHL) + 12.0 * (trH2 * trH2 + trL2 * trL2) - 8.0 * (trH2 * trL2);

        arma::cx_double action = trD4 + this->sim_data.g2 * trD2;
        if (action.imag() > 1e-9)
        {
            throw std::runtime_error("Action has an imaginary component greater than tolerance.");
        }

        return action.real();
    }
};

#endif // TYPE11SIMULATION_H