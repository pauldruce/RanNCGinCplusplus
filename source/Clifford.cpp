#include "Clifford.h"
using namespace std;
using namespace arma;

Clifford::Clifford(int p, int q)
{
   this->p = p;
   this->q = q;
   
   type = std::make_pair(p,q);
   n = p + q;
   k = (n%2==0) ? n/2 : (n-1)/2;
   s = (q - p) % 8;
   ssp1 = s * (s + 1);
   matrix_size = pow(2, k);
   setup();
}

/* get_chirality():
 
 Once we have constructed all of the clifford generators, they are stored in
 "generators" and are ordered so the hermitian generators first and the
 anti-hermitian generators are second.
 
 We can then calculate the Chirality operator using: chirality = i^{s(s+1)/2}
 gamma^1 gamma^2 ... gamma^n
 */
cx_mat Clifford::get_chirality()
{
   cx_mat c;
   c = eye<cx_mat>(matrix_size, matrix_size) ;
   c *= std::pow(std::complex<double>(0, 1), ssp1/2);
   for (auto gam : generators)
      c *= gam;
   c.clean(datum::eps);
   return c;
   
}

cx_mat Clifford::generate_chirality(vector<cx_mat>& input_generators)
{
   cx_mat c = eye<cx_mat>(matrix_size,matrix_size);
   c *= std::pow(std::complex<double>(0, 1), ssp1 / 2);
   for (auto gam : input_generators)
      c *= gam;
   c.clean(datum::eps);
   return c;
   
}

/*  introduce():
 
 This function prints out the type, the expected matrix dimension, the K0
 dimension, the generators and the chirality operator. This is useful to
 check that everything is working as planned.
 TODO: Check that the generators satisfy the right anti-commutator
 relationships.
 
 */

void Clifford::introduce()
{
   cout<< "My type is: ("<< p <<","<< q << ")\n";
   cout<<"My matrices will be"<< matrix_size<< "x" <<matrix_size<< ", and I have K0 dimension s ="<< s <<endl;
      // Print out the generators
   cout << "The generators are:" << endl;
   int i = 0;
   for( auto gen: generators)
   {
      // Starts the gamma labels from one.
      cout << "Gamma_"<<i<< "=" << endl;
      gen.print();
      cout << endl;
      i++;
   }
   cout<< "The chirality operator is:"<< endl;
   chirality.print();
}

/* setup():
 
 This is the main function of the class which setups up the generators.
 The small dimensional cases of (1,0), (0,1), (1,1), (2,0) and (0,2) are
 coded in by hand. The higher dimensional cases are then constructed using
 the product mechanism that can be found in Lawson and Michelsohn for
 instance, or there is a more to the point demonstration in my PhD thesis.
 
 Note that for type (1,3) this does not produce the gamma matrices in the
 Dirac basis or the Chiral/Weyl basis or the Majorana basis.
 
 TODO: There is a function to apply a function that converts it to the
 Dirac/Chiral/Majorana basis
 */

void Clifford::setup()
{
   // Type (0,0) setup
   if (p == 0 and q == 0)
   {
      generators.clear();
      chirality = (cx_mat){0};
   }
   // Type (0,1) setup
   else if (p == 0 and q == 1)
   {
      cx_mat gamma1 = {std::complex<double>(0, 1)};
      generators = {gamma1};
      chirality = generate_chirality(generators);
   }
   // Type (1,0) setup
   else if (p == 1 and q == 0)
   {
      cx_mat gamma1 = {std::complex<double>(1, 0)};
      generators = {gamma1};
      chirality = generate_chirality(generators);
   }
   // Type (0,2) setup
   else if (p == 0 and q == 2)
   {
      cx_mat gamma1 = {{complex<double>(0, 1), 0}, {0, complex<double>(0, -1)}};
      cx_mat gamma2 = {{0, complex<double>(1, 0)}, {complex<double>(-1, 0), 0}};
      generators = {gamma1, gamma2};
      chirality = generate_chirality(generators);
   }
   // Type (1,1) setup
   else if (p == 1 and q == 1)
   {
      cx_mat gamma1 = {{complex<double>(1, 0), 0}, {0, complex<double>(-1, 0)}};
      cx_mat gamma2 = {{0, complex<double>(1, 0)}, {complex<double>(-1, 0), 0}};
      generators = {gamma1, gamma2};
      chirality = generate_chirality(generators);
   }
   // Type (2,0) setup
   else if (p == 2 and q == 0)
   {
      cx_mat gamma1 = {{complex<double>(1, 0), 0}, {0, complex<double>(-1, 0)}};
      cx_mat gamma2 = {{0, complex<double>(1, 0)}, {complex<double>(1, 0), 0}};
      generators = {gamma1, gamma2};
      chirality = generate_chirality(generators);
   }
   // This is the code for all Clifford Algebras/Modules with p+q>2. This
   // procedure is outlined in Lawson, H.B., Jr, Michelsohn, M.-L.: Spin
   // Geometry. Princeton University Press, Princeton, NJ (1989) and a more
   // direct and to the point description of it is in my thesis which can be
   // found at my website: https://pauldruce.github.io/docs/
   if (n > 2)
   {
      int d_p = (int)p / 2; // Number of times to product with (2,0)
      int d_q = (int)q / 2; // Number of times to product with (0,2)
                            // Number of times to product with (1,0) i.e. either once, or not at all
      int r_p = p % 2;
      // Number of times to product with (0,1) i.e. either once, or not at all
      int r_q = q % 2;
      
      /*
       There are four options:
       -  p = even and q = even
       -  p = even and q = odd
       -  p = odd and q = even
       -  p = odd and q = odd
       
       If p = even and q = even, Then start with (0,2) or (2,0) and we can
       product (2,0) and (0,2) together the correct number of times.
       
       If p = even and q = odd, start with (0,1) and product with (2,0) or (0,2)
       the correct number of times
       
       If p = odd and q = even, start with (1,0) and product with (2,0) or (0,2)
       the correct number of times
       
       If p = odd and q = odd, start with (1,1) and product with (2,0) or (0,2)
       
       */
      
      auto cliff02 = Clifford(0, 2);
      auto cliff20 = Clifford(2, 0);
      auto cliff11 = Clifford(1, 1);
      vector<cx_mat> input_gens;
      vector<cx_mat> holder_gens;
      
      // If p = even and q = even
      if (r_p == 0 and r_q == 0)
      {
         if (d_p != 0 and d_q == 0)
         {
            input_gens = cliff20.generators;
            d_p = d_p - 1; // Adjust for the fact we start from (2,0)
            module_dim = 2;
         }
         if (d_p == 0 and d_q != 0)
         {
            input_gens = cliff02.generators;
            d_q = d_q - 1; // Adjust for the fact we start from (0,2)
            module_dim = 2;
         }
         if (d_p != 0 and d_q != 0)
         { // This branch needs thought
            input_gens = cliff20.generators;
            d_p = d_p - 1; // Adjust for the fact we start from (2,0)
            module_dim = 2;
         }
      }
      // If p = even and q = odd
      else if (r_p == 0 and r_q == 1)
      {
         input_gens = Clifford(0, 1).generators;
         module_dim = 1;
      }
      // If p = odd and q = even
      else if (r_p == 1 and r_q == 0)
      {
         input_gens = Clifford(1, 0).generators;
         module_dim = 1;
      }
      // If p = odd and q = odd
      else if (r_p == 1 and r_q == 1)
      {
         input_gens = Clifford(1, 1).generators;
         module_dim = 2;
      }
      /*
       Products with (2,0) and (0,2)
       
       This procedure requires the first Clifford module to have even s. As both
       Cliff(2,0) and Cliff(0,2) have even s. This works.
       
       We will do the following
       let M = M(0,0)
       do M = M(2,0) x M  d_p times
       do M = M(0,2) x M  d_q times
       so that M = M(2*d_p, 2*d_q)
       */
      
      for (int i = 0; i < d_p; i++)
      {
         
         // This gives me the matrix dimension of the second module we are
         // producting.
         for (auto gen : cliff20.generators)
         {
            holder_gens.push_back(kron(gen, eye<cx_mat>(module_dim, module_dim)));
         }
         for (auto gen : input_gens)
         {
            holder_gens.push_back(kron(cliff20.chirality, gen));
            
         }
         // After producting with a (2,0) the module_dim =
         // matrix_dimension is increased by a factor of 2
         module_dim *= 2;
         // Set the input_gens to be the just constructed
         // generators so we can then repeat the process if
         // necessary.
         input_gens = holder_gens;
         
      }
      
      for (int i = 0; i < d_q; i++)
      {
         // This gives me the matrix dimension of the second module we are producting.
         for (auto gen : cliff02.generators)
         {
            holder_gens.push_back(kron(gen, eye<cx_mat>(module_dim, module_dim)));
         }
         for (auto gen : input_gens)
         {
            
            holder_gens.push_back(kron(cliff02.chirality, gen));
         }
         // After producting with a (2,0) the module_dim = matrix_dimension is increased by a factor of 2
         module_dim *= 2;
         // Set the input_gens to be the just constructed generators so we can then repeat the process if necessary.
         input_gens = holder_gens;
      }
      
      // Order the generators so the Hermitian elements are first, then the anti-Heritian
      vector<cx_mat> herm;
      vector<cx_mat> anti_herm;
      for (auto gen : input_gens)
      {
         gen.clean(datum::eps);
         if (gen.is_hermitian()) // np.array_equal(gen.conj().T, gen) == True
         {
            herm.push_back(gen);
         }
         else
         {
            anti_herm.push_back(gen);
         }
      }
      
      // Set the generators to what has been calculated.
      // The + operation here is action of Python lists, so it concatenates
      // them. It doesn't add the operators together.
      generators = herm;
      generators.insert(generators.end(), anti_herm.begin(), anti_herm.end());
      //      generators.insert(generators.end(), herm.begin(), herm.end());
      //      generators.insert(generators.end(), anti_herm.begin(), anti_herm.end());
      chirality = get_chirality();
   }
}





