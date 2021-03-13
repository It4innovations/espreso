#ifndef SRC_WRAPPERS_CSPARSE_W_CSPARSE_H_
#define SRC_WRAPPERS_CSPARSE_W_CSPARSE_H_

#include "feti/generic/SparseMatrix.h"

namespace espreso {
namespace csparse {

// DEPRECATED - Will be removed
// CSparse factorization - 0-based indexing, CSC-format
// Factorize A on CPU, solve mrhs linear system and SpMM on GPU
// void CreateLscGpu(SparseMatrix& A, SparseMatrix& B, int order, int tol, int gpu_id, int print_output, SparseMatrix& SC);

// Performs Cholesky factorization and keeps factors
// Important: factors must be freed via FreeCholFactor function
// A is square/symmetric CSR matrix in 1-based indexing
// order 0 = natural, 1 = amd(A+A')
// L is output triangular matrix in CSC
// perm is output permutation vector
void FactorizeChol(SparseMatrix& A, int order, int& L_nnz, int*& L_row_indexes, int*& L_col_pointers, double*& L_values, int*& perm);

void FreeCholFactor(int* L_row_indexes, int* L_col_pointers, double* L_values, int* perm);

// Performs LU factorization and keeps factors
// Important: factors must be freed via FreeLuFactors function
// A is square/symmetric CSR matrix in 1-based indexing
// order 0 = natural, 1 = amd(A+A'), 2 = amd(S'*S), 3 = amd(A'*A)
// L is output lower triangular matrix in CSC
// U is output upper triangular matrix in CSC
// perm is output permutation vector
// perm_2 is output permutation vector
void FactorizeLu(SparseMatrix& A, int order, int& L_nnz, int*& L_row_indexes, int*& L_col_pointers, double*& L_values,
 int& U_nnz, int*& U_row_indexes, int*& U_col_pointers, double*& U_values, int*& perm, int*& perm_2);

void FreeLuFactors(int* L_row_indexes, int* L_col_pointers, double* L_values,
 int* U_row_indexes, int* U_col_pointers, double* U_values, int* perm, int* perm_2);
 
}
}

#endif /* SRC_WRAPPERS_CSPARSE_W_CSPARSE_H_ */
