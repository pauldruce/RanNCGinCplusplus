#include <doctest/doctest.h>
#include "Type11Simulation.hpp"

TEST_SUITE("Type11Simulation Class Tests - Testing overridden action method of Simulation class")
{
    TEST_CASE("Type11Simulation Class Tests - Testing overridden action method of Simulation class")
    {
        SUBCASE("Throws error if incorrect type Dirac operator is passed")
        {
            int N = 5;
            Clifford type12(1, 2);
            DiracOperator dirac(type12, N);
            SimulationData sim_data;

            CHECK_THROWS_AS(Type11Simulation sim(dirac, sim_data), std::runtime_error);
        }
    }
}