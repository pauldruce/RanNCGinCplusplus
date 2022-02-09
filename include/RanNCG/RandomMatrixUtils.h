//
//  RandomMatrixUtils.h
//  RandomNCG
//
//  Created by Paul Druce on 29/12/2021.
//

#ifndef RandomMatrixUtils_h
#define RandomMatrixUtils_h

#include <armadillo>

// Create a NxN complex matrix with entires uniformly selected from [-1,1] + i[-1,1]
arma::cx_mat random_complex(int& N);

// Creates an NxN element of a Gaussian Hermitian Ensemble
arma::cx_mat random_hermitian(int& N);


// Returns the commutator of a matrix as [M,-] = M x I - I x M.T
arma::cx_mat comm(arma::cx_mat& M);

// Returns the anti-commutator of an n-by-n complex matrix: {M,-} = M \otimes I + I\otimes M.T
arma::cx_mat anticomm(arma::cx_mat& M);

#endif /* RandomMatrixUtils_h */
