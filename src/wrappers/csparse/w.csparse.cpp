#include "w.csparse.h"

#ifdef HAVE_CSPARSE
#include "cs.h"
#include "esinfo/eslog.h"

namespace espreso {

void csparse::CreateLscGpu(SparseMatrix& A, SparseMatrix& B, int order, int tol, int gpu_id, int print_output, SparseMatrix& SC) {
    // | A (B)t |
    // | B   0  |

    int fact = 0; // 0 = LU, 1 = Cholesky
    if(A.type == 'S') {
        fact = 1;
    }

    // Convert 1-based to 0-based indexing
    for(int& i : A.CSR_I_row_indices) {
        i--;
    }
    for(int& i : A.CSR_J_col_indices) {
        i--;
    }

    cs* A_cs = (cs*) cs_calloc (1, sizeof (cs));
	// A_cs in CSC format
    A_cs->m = A.rows;
    A_cs->n = A.cols;
    A_cs->p = A.CSR_I_row_indices.data();
    A_cs->i = A.CSR_J_col_indices.data();
    A_cs->x = A.CSR_V_values.data();
    A_cs->nz = -1; // # of entries in triplet matrix, -1 for compressed-col
    A_cs->nzmax = A.CSR_V_values.size();

    // Convert 1-based to 0-based indexing
    for(int& i : B.CSR_I_row_indices) {
        i--;
    }
    for(int& i : B.CSR_J_col_indices) {
        i--;
    }

    cs* B_cs = (cs*) cs_calloc (1, sizeof (cs)) ;
    // B_cs in CSR format
    B_cs->m = B.rows;
    B_cs->n = B.cols;
    B_cs->p = B.CSR_I_row_indices.data();
    B_cs->i = B.CSR_J_col_indices.data();
    B_cs->x = B.CSR_V_values.data();
    B_cs->nz = -1; // # of entries in triplet matrix, -1 for compressed-col
    B_cs->nzmax = B.CSR_V_values.size();

	// Factorize K on CPU, solve mrhs linear system and SpMM on GPU
	int ok = create_lsc(fact, order, A_cs, B_cs, tol, SC.d_dense_values, gpu_id, print_output);

    if(!ok) eslog::error("ESPRESO run-time error: CSparse library factorization failed.\n");

    SC.cols = B.rows;
    SC.rows = B.rows;
    SC.type = 'G';

    cs_free(B_cs);
    cs_free(A_cs);
}

}

#endif