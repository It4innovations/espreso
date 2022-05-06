#include "w.csparse.h"
#include "esinfo/eslog.h"

//#ifndef HAVE_CSPARSE
//
//namespace espreso {
//
//// void csparse::CreateLscGpu(SparseMatrix& A, SparseMatrix& B, int order, int tol, int gpu_id, int print_output, SparseMatrix& SC) {
////     eslog::warning("ESPRESO run-time error: cannot call CSparse library (the library is not linked).\n");
//// }
//
//void csparse::FactorizeChol(SparseMatrix& A, int order, int& L_nnz, int*& L_row_indexes, int*& L_col_pointers, double*& L_values, int*& perm) {
//    eslog::warning("ESPRESO run-time error: cannot call CSparse library (the library is not linked).\n");
//}
//
//void csparse::FreeCholFactor(int* L_row_indexes, int* L_col_pointers, double* L_values, int* perm) {
//    eslog::warning("ESPRESO run-time error: cannot call CSparse library (the library is not linked).\n");
//}
//
//void csparse::FactorizeLu(SparseMatrix& A, int order, int& L_nnz, int*& L_row_indexes, int*& L_col_pointers, double*& L_values,
// int& U_nnz, int*& U_row_indexes, int*& U_col_pointers, double*& U_values, int*& perm, int*& perm_2) {
//    eslog::warning("ESPRESO run-time error: cannot call CSparse library (the library is not linked).\n");
//}
//
//void csparse::FreeLuFactors(int* L_row_indexes, int* L_col_pointers, double* L_values,
// int* U_row_indexes, int* U_col_pointers, double* U_values, int* perm, int* perm_2) {
//    eslog::warning("ESPRESO run-time error: cannot call CSparse library (the library is not linked).\n");
//}
//
//}
//
//#endif
