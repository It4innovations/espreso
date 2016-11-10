/*! \file   DissectionDefault.hpp
    \brief  definition of default value for factorization
    \author Atsushi Suzuki, Laboratoire Jacques-Louis Lions
    \date   Jul. 24th 2015
    \date   Sep. 29th 2015
    \date   Feb. 29th 2016
*/

#ifndef _DISSECTION_DEFAULT_
#define _DISSECTION_DEFAULT_

#define SCOTCH_DECOMPOSER  0
#define METIS_DECOMPOSER   1
#define TRIDIAG_DECOMPOSER 2

#define NO_SCALING       0  // needs to be compatible to the definition in
#define DIAGONAL_SCALING 1  // SparseMatrix<T, U>::normalize(), SparseMatrix.cpp
#define KKT_SCALING      2  //

#define MINNODES     256    // minimum size of the first layer of dissection
#define SIZE_TRIDIAG 1000   // more than this value, dissection is used

#define DIM_AUG_KERN     4        // appropriate for indefinite matrix 
#define EPS_PIVOT        1.0e-2

#define TOL_PIVOT        1.0e-5   // for recursion of sparse factorization
#define MIN_TRIDIAG_SIZE 50       // the size to avoid Cuthill-McKee ordering

#endif
