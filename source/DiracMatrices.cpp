#include "DiracMatrices.h"
#include <armadillo>

using namespace std;
using namespace arma;

DiracMatrices::DiracMatrices()
    : Clifford(1, 3)
{

   // Define the dirac basis for type (1,3)
   gamma0 = {{0, 0, 1, 0},
             {0, 0, 0, 1},
             {1, 0, 0, 0},
             {0, 1, 0, 0}};
   gamma1 = {{0, 0, 0, 1},
             {0, 0, 1, 0},
             {0, -1, 0, 0},
             {-1, 0, 0, 0}};
   gamma2 = {{0, 0, 0, complex<double>(0, -1)},
             {0, 0, complex<double>(0, 1), 0},
             {0, complex<double>(0, 1), 0, 0},
             {complex<double>(0, -1), 0, 0, 0}};
   gamma3 = {{0, 0, 1, 0},
             {0, 0, 0, -1},
             {-1, 0, 0, 0},
             {0, 1, 0, 0}};

   generators = {gamma0, gamma1, gamma2, gamma3};
   dirac_gens = generators;
   dirac_chirality = generate_chirality(dirac_gens);

   // Define the majorana basis
   cx_mat majorana0 = {{0, 0, 0, complex<double>(0, -1)},
                       {0, 0, complex<double>(0, 1), 0},
                       {0, complex<double>(0, -1), 0, 0},
                       {complex<double>(0, 1), 0, 0, 0}};
   cx_mat majorana1 = {{complex<double>(0, 1), 0, 0, 0},
                       {0, complex<double>(0, -1), 0, 0},
                       {0, 0, complex<double>(0, 1), 0},
                       {0, 0, 0, complex<double>(0, -1)}};
   cx_mat majorana2 = {{0, 0, 0, complex<double>(0, 1)},
                       {0, 0, complex<double>(0, -1), 0},
                       {0, complex<double>(0, -1), 0, 0},
                       {complex<double>(0, 1), 0, 0, 0}};
   cx_mat majorana3 = {{0, complex<double>(0, -1), 0, 0},
                       {complex<double>(0, -1), 0, 0, 0},
                       {0, 0, 0, complex<double>(0, -1)},
                       {0, 0, complex<double>(0, 1), 0}};
   majorana_gens = {majorana0, majorana1, majorana2, majorana3};
   majorana_chirality = -1 * generate_chirality(majorana_gens);

   // Define the Chiral basis
   cx_mat chiral0 = {{0, 0, 1, 0},
                     {0, 0, 0, 1},
                     {1, 0, 0, 0},
                     {0, 1, 0, 0}};
   cx_mat chiral1 = {{0, 0, 0, 1},
                     {0, 0, 1, 0},
                     {0, 0, 0, -1},
                     {0, 0, -1, 0}};
   cx_mat chiral2 = {{0, 0, 0, complex<double>(0, -1)},
                     {0, 0, complex<double>(0, 1), 0},
                     {0, complex<double>(0, 1), 0, 0},
                     {complex<double>(0, -1), 0, 0, 0}};
   cx_mat chiral3 = {{0, 0, 1, 0}, {0, 0, 0, -1}, {-1, 0, 0, 0}, {0, 1, 0, 0}};

   chiral_gens = {chiral0, chiral1, chiral2, chiral3};
   chiral_chirality = generate_chirality(chiral_gens);

   // Get the new chirality operator
   chirality = get_chirality();
}

void DiracMatrices::introduce()
{
   cout << "My type is: (" << p << "," << q << ")\n";
   cout << "My matrices will be " << matrix_size << "x" << matrix_size << ", and I have K0 dimension s=" << s << endl;

   // Print out the generators
   cout << "The generators are:";
   int i = 0;
   for (auto gen : generators)
   {
      // Starts the gamma label from zero
      cout << "Gamma_" << i << "= " << endl;
      gen.print();
      cout << endl;
      i++;
   }

   cout << "The chirality operator is:" << endl;
   chirality.print();
}
