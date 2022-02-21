//
//  RandomMatrixUtils.cpp
//  RandomNCG
//
//  Created by Paul Druce on 29/12/2021.
//

#include "RandomMatrixUtils.h"
//#include <cassert>
using namespace arma;
using namespace std;

cx_mat comm(cx_mat& M)
{
   /*
    * returns the commutator of a matrix as [M,-] = M x I - I x M.T
    *
    *  To express the right action of a Matrix onto another, we often write M_n(C) as C^n \otimes (C^n)^*,
    *  where the * represents the dual vector space. In terms of vectors, if we view C^n as column vectors,
    *  then (C^n)^* are row vectors, so v^* = v^T. So a matrix acts on the left of an element v \otimes v^*, by Mv \otimes I.v^*,
    *  where I is the identity matrix. So a matrix acts on the right by I.v \otimes M^T.v^*
    */

   if (!M.is_square())
   {
      cout << "Trying to form commutator with non-square matrix!" << endl;
      cout << "Returning input matrix unchanged." << endl;
      throw "Input matrix is not square and you can not form commutator with a non-square matrix.";
   }

   cx_mat I = eye<cx_mat>(size(M));
   return kron(M, I) - kron(I, M.st());

}


cx_mat anticomm(cx_mat& M)
{
   /*
    * Returns the anti-commutator of an n-by-n complex matrix: {M,-} = M \otimes I + I\otimes M.T
    * To express the right action of a Matrix onto another, we often write M_n(C) as C^n \otimes (C^n)^*,
    * where the * represents the dual vector space. In terms of vectors, if we view C^n as column vectors,
    * then (C^n)^* are row vectors, so v^* = v^T. So a matrix acts on the left of an element v \otimes v^*, by Mv \otimes I.v^*,
    * where I is the identity matrix. So a matrix acts on the right by I.v \otimes M^T.v^*
    *
    * So this function returns {M,-} = M \otimes I + I \otimes M.T
    */
   
   if (!M.is_square())
   {
      cout << "Trying to form commutator with non-square matrix!" << endl;
      cout << "Returning input matrix unchanged." << endl;
      throw "Input matrix is not square and you can not form anti-commutator with a non-square matrix.";
   }
   cx_mat I = eye<cx_mat>(size(M));
   return kron(M, I) + kron(I, M.st());
}

cx_mat random_complex(int& N)
{
   mat A(N,N,fill::randu);
   mat B(N,N,fill::randu);
   
   // Need elements between [-1,1] and randu gives values between [0,1]
   A = 2*A - 1;
   B = 2*B - 1;
   return cx_mat(A, B);
}

cx_mat random_hermitian(int& N)
{  
   mat A(N, N, fill::randn);
   mat B(N, N, fill::randn);
   cx_mat C = cx_mat(A, B); // Creates a matrix A + iB
   
   C = C + C.t(); // In Armadillo, t() is the Hermitian transpose (the adjoint) and st() is the ordinary transpose for complex matrices.

   assert(C.is_hermitian());

   return C;
}
