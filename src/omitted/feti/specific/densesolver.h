
#ifndef SOLVER_SPECIFIC_DENSESOLVER_H_
#define SOLVER_SPECIFIC_DENSESOLVER_H_

// #include "esbasis.h"
// #include "generic/utils.h"
#include "feti/generic/SparseMatrix.h"

namespace espreso {

class DenseSolver
{

public:
    virtual ~DenseSolver() {};

    virtual void ImportMatrix    (SparseMatrix & A) = 0;
    virtual void ImportMatrix_fl (SparseMatrix & A) = 0;

    virtual void ImportMatrix_wo_Copy (SparseMatrix & A) = 0;
    virtual void ImportMatrix_wo_Copy_fl(SparseMatrix & A) = 0;

    virtual int Factorization(const std::string &str) = 0;
    virtual void Clear() = 0;
    virtual void SetThreaded() = 0;


    virtual void Solve( SEQ_VECTOR <double> & rhs, SEQ_VECTOR <double> & sol, MKL_INT rhs_start_index, MKL_INT sol_start_index) = 0;
    virtual void Solve( SEQ_VECTOR <double> & rhs, SEQ_VECTOR <double> & sol, MKL_INT n_rhs) = 0;
    virtual void Solve( SEQ_VECTOR <double> & rhs_sol) = 0;

    virtual void SolveMat_Sparse( SparseMatrix & A ) = 0;
    virtual void SolveMat_Sparse( SparseMatrix & A_in, SparseMatrix & B_out ) = 0;
    virtual void SolveMat_Sparse( SparseMatrix & A_in, SparseMatrix & B_out, char T_for_input_matrix_is_transposed_N_input_matrix_is_NOT_transposed ) = 0;

    virtual void SolveMat_Dense( SparseMatrix & A_in_out ) = 0;
    virtual void SolveMat_Dense( SparseMatrix & A_in, SparseMatrix & B_out ) = 0;

    virtual void SolveCG(SparseMatrix & A_in, SEQ_VECTOR <double> & rhs_in, SEQ_VECTOR <double> & sol, SEQ_VECTOR <double> & initial_guess) = 0;
    virtual void SolveCG(SparseMatrix & A_in, std::vector <double> & rhs, std::vector <double> & sol) = 0;
    virtual void SolveCG(SparseMatrix & A_in, std::vector <double> & rhs_sol) = 0;

};

}

#endif
