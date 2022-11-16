#include "w.csparse.h"

#ifdef HAVE_CSPARSE
#include "cs.h"
#include "esinfo/eslog.h"

namespace espreso {

 /* 1 if A is square & upper tri., -1 if square & lower tri., 0 otherwise */
int IsSym(cs *A) {
    int is_upper, is_lower, j, p, n = A->n, m = A->m, *Ap = A->p, *Ai = A->i ;
    if (m != n) return (0) ;
    is_upper = 1 ;
    is_lower = 1 ;
    for (j = 0 ; j < n ; j++)
    {
        for (p = Ap [j] ; p < Ap [j+1] ; p++)
        {
            if (Ai [p] > j) is_upper = 0 ;
            if (Ai [p] < j) is_lower = 0 ;
        }
    }
    return (is_upper ? 1 : (is_lower ? -1 : 0)) ;
}

/* true for off-diagonal entries */
int Dropdiag (int i, int j, double aij, void *other) {
    return (i != j) ;
}


/* C = A + triu(A,1)' */
cs *MakeSym (cs *A) {
    cs *AT, *C ;
    AT = cs_transpose (A, 1) ;          /* AT = A' */
    cs_fkeep (AT, &Dropdiag, NULL) ;    /* drop diagonal entries from AT */
    C = cs_add (A, AT, 1, 1) ;          /* C = A+AT */
    cs_spfree (AT) ;
    return (C) ;
}

// DEPRECATED - to be removed
// void csparse::CreateLscGpu(SparseMatrix& A, SparseMatrix& B, int order, int tol, int gpu_id, int print_output, SparseMatrix& SC) {
//     // | A (B)t |
//     // | B   0  |

//     int fact = 0; // 0 = LU, 1 = Cholesky
//     if(A.type == 'S') {
//         fact = 1;
//         eslog::error("ESPRESO run-time error: CSparse Chol factorization wrapper - triangular matrix handling not implemented yet.\n");
//     }

//     // Convert 1-based to 0-based indexing
//     for(int& i : A.CSR_I_row_indices) {
//         i--;
//     }
//     for(int& i : A.CSR_J_col_indices) {
//         i--;
//     }

//     cs* A_cs = (cs*) cs_calloc (1, sizeof (cs));
// 	// A_cs in CSC format
//     A_cs->m = A.rows;
//     A_cs->n = A.cols;
//     A_cs->p = A.CSR_I_row_indices.data();
//     A_cs->i = A.CSR_J_col_indices.data();
//     A_cs->x = A.CSR_V_values.data();
//     A_cs->nz = -1; // # of entries in triplet matrix, -1 for compressed-col
//     A_cs->nzmax = A.CSR_V_values.size();

//     // Convert 1-based to 0-based indexing
//     for(int& i : B.CSR_I_row_indices) {
//         i--;
//     }
//     for(int& i : B.CSR_J_col_indices) {
//         i--;
//     }

//     cs* B_cs = (cs*) cs_calloc (1, sizeof (cs)) ;
//     // B_cs in CSR format
//     B_cs->m = B.rows;
//     B_cs->n = B.cols;
//     B_cs->p = B.CSR_I_row_indices.data();
//     B_cs->i = B.CSR_J_col_indices.data();
//     B_cs->x = B.CSR_V_values.data();
//     B_cs->nz = -1; // # of entries in triplet matrix, -1 for compressed-col
//     B_cs->nzmax = B.CSR_V_values.size();

// 	// Factorize K on CPU, solve mrhs linear system and SpMM on GPU
// 	int ok = create_lsc(fact, order, A_cs, B_cs, tol, SC.d_dense_values, gpu_id, print_output);

//     if(!ok) eslog::error("ESPRESO run-time error: CSparse library factorization failed.\n");

//     SC.cols = B.rows;
//     SC.rows = B.rows;
//     SC.type = 'G';

//     cs_free(B_cs);
//     cs_free(A_cs);
// }


void csparse::FactorizeChol(SparseMatrix& A, int order, int& L_nnz, int*& L_row_indexes, int*& L_col_pointers, double*& L_values, int*& perm) {
    // Convert 1-based to 0-based indexing
    for(int& i : A.CSR_I_row_indices) {
        i--;
    }
    for(int& i : A.CSR_J_col_indices) {
        i--;
    }
    // Import SparseMatrix
    cs* A_cs = (cs*) cs_calloc(1, sizeof(cs));
    A_cs->m = A.rows;
    A_cs->n = A.cols;
    A_cs->p = A.CSR_I_row_indices.data();
    A_cs->i = A.CSR_J_col_indices.data();
    A_cs->x = A.CSR_V_values.data();
    A_cs->nz = -1; // # of entries in triplet matrix, -1 for compressed-col
    A_cs->nzmax = A.CSR_V_values.size();

    int sym = IsSym(A_cs);
    if (sym == -1)
        A_cs = cs_transpose (A_cs, 1);

    // C = A + triu(A,1)', or C=A
//    cs* A_cs_sym = sym ? MakeSym(A_cs) : A_cs;

    // Ordering and symbolic analysis
    css *S = cs_schol(order, A_cs);

    // Numeric Cholesky factorization
    csn *N = cs_chol(A_cs, S);

    // Get L factor
    L_nnz = N->L->nzmax;
    L_row_indexes = N->L->i;
    L_col_pointers = N->L->p;
    L_values = N->L->x;

    // Get ordering (permutation vector)
    perm = S->pinv ? S->pinv : NULL;

    // Convert 0-based to 1-based indexing
    for(int& i : A.CSR_I_row_indices) {
        i++;
    }
    for(int& i : A.CSR_J_col_indices) {
        i++;
    }

    // Prevent deleting data
    S->pinv = NULL;
    N->L->i = NULL;
    N->L->p = NULL;
    N->L->x = NULL;

    // Free memory
    cs_free(A_cs);
    cs_nfree(N);
    cs_sfree(S);
}

void csparse::FreeCholFactor(int* L_row_indexes, int* L_col_pointers, double* L_values, int* perm) {
    cs_free(L_row_indexes);
    cs_free(L_col_pointers);
    cs_free(L_values);
    cs_free(perm);
}

void csparse::FactorizeLu(SparseMatrix& A, int order, int& L_nnz, int*& L_row_indexes, int*& L_col_pointers, double*& L_values,
 int& U_nnz, int*& U_row_indexes, int*& U_col_pointers, double*& U_values, int*& perm, int*& perm_2) {
    // A_cs in CSC format so transpose A in CSR
    A.MatTranspose();
    
    // Convert 1-based to 0-based indexing
    for(int& i : A.CSR_I_row_indices) {
        i--;
    }
    for(int& i : A.CSR_J_col_indices) {
        i--;
    }

    // Import SparseMatrix
    cs *A_cs = (cs*) cs_calloc(1, sizeof(cs));
    A_cs->m = A.rows;
    A_cs->n = A.cols;
    A_cs->p = A.CSR_I_row_indices.data();
    A_cs->i = A.CSR_J_col_indices.data();
    A_cs->x = A.CSR_V_values.data();
    A_cs->nz = -1; // # of entries in triplet matrix, -1 for compressed-col
    A_cs->nzmax = A.CSR_V_values.size();

    int tol = 1;
    // TODO: if sym -> should be possible to symmetrize like in Factorize Chol in case of LU for symmetric matrix
    // if(sym) {
    //     tol = 0.001;
    // }

    // Ordering and symbolic analysis
    css *S = cs_sqr(order, A_cs, 0);

    // Numeric Cholesky factorization
    csn *N = cs_lu(A_cs, S, tol);

    // Get L and U factors
    L_nnz = N->L->nzmax;
    L_row_indexes = N->L->i;
    L_col_pointers = N->L->p;
    L_values = N->L->x;

    U_nnz = N->U->nzmax;
    U_row_indexes = N->U->i;
    U_col_pointers = N->U->p;
    U_values = N->U->x;

    // Get ordering (permutation vector)
    perm = N->pinv ? N->pinv : NULL;
    perm_2 = S->q ? S->q : NULL;

    // Convert 0-based to 1-based indexing
    for(int& i : A.CSR_I_row_indices) {
        i++;
    }
    for(int& i : A.CSR_J_col_indices) {
        i++;
    }

    // Convert back A from CSC to CSR
    A.MatTranspose();

    // Prevent deleting data
    N->pinv = NULL;
    S->q = NULL;

    N->L->i = NULL;
    N->L->p = NULL;
    N->L->x = NULL;

    N->U->i = NULL;
    N->U->p = NULL;
    N->U->x = NULL;

    // Free memory
    cs_free(A_cs);
    cs_nfree(N);
    cs_sfree(S);
}

void csparse::FreeLuFactors(int* L_row_indexes, int* L_col_pointers, double* L_values,
 int* U_row_indexes, int* U_col_pointers, double* U_values, int* perm, int* perm_2) {
    cs_free(L_row_indexes);
    cs_free(L_col_pointers);
    cs_free(L_values);
    cs_free(U_row_indexes);
    cs_free(U_col_pointers);
    cs_free(U_values);
    cs_free(perm);
    cs_free(perm_2);
}

}

#endif
