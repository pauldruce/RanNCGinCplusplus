#pragma once
#include <armadillo>
#include <utility>

class Clifford {
public:
   Clifford(){};
   Clifford(int p, int q);
protected:
   // Define constants for the Clifford type (p,q), the size n=p+4, and the matrx
   // dimensions 2^k, where k  = n/2 or (n-1)/2 depending on n. Also have s(s+1)
   // defined where s = q-p
   int p, q,n, k, matrix_size, s, ssp1;
   int module_dim;
   void setup();
   arma::cx_mat generate_chirality(std::vector<arma::cx_mat>& input_generators);
   
public:
   std::vector<arma::cx_mat> generators;
   arma::cx_mat chirality;
   std::pair<int,int> type;
   virtual void introduce();
   arma::cx_mat get_chirality();

   inline int get_matrix_size() { return this->matrix_size; };
   inline int get_number_of_generators(){return n;};
   inline int get_p(){return p;};
   inline int get_q(){return q;};
};
