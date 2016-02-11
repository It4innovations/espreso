
#ifndef SOLVER_SPECIFIC_SPARSESOLVER_H_
#define SOLVER_SPECIFIC_SPARSESOLVER_H_

#include "../generic/utils.h"
#include "../generic/SparseMatrix.h"

class SparseSolver
{

public:
	virtual ~SparseSolver() {};

	virtual void ImportMatrix(SparseMatrix & A) = 0;
	virtual void ImportMatrix_fl(SparseMatrix & A) = 0;

	virtual void ImportMatrix_wo_Copy(SparseMatrix & A) = 0;

	virtual void Factorization(const std::string &str) = 0;
	virtual void Clear() = 0;
	virtual void SetThreaded() = 0;


	virtual void Solve( SEQ_VECTOR <double> & rhs, SEQ_VECTOR <double> & sol, MKL_INT rhs_start_index, MKL_INT sol_start_index) = 0;
	virtual void Solve( SEQ_VECTOR <double> & rhs, SEQ_VECTOR <double> & sol, MKL_INT n_rhs) = 0;
	virtual void Solve( SEQ_VECTOR <double> & rhs_sol) = 0;

	virtual void SolveMat_Dense( SparseMatrix & A_in_out ) = 0;
	virtual void SolveMat_Dense( SparseMatrix & A_in, SparseMatrix & B_out ) = 0;

	virtual void SolveMatF( SparseMatrix & A_in, SparseMatrix & B_out, bool isThreaded ) = 0;

	virtual void SolveMat_Sparse( SparseMatrix & A ) = 0;
	virtual void SolveMat_Sparse( SparseMatrix & A_in, SparseMatrix & B_out ) = 0;
	virtual void SolveMat_Sparse( SparseMatrix & A_in, SparseMatrix & B_out, char T_for_input_matrix_is_transposed_N_input_matrix_is_NOT_transposed ) = 0;

	virtual void Create_SC( SparseMatrix & B_out, MKL_INT sc_size, bool isThreaded ) = 0;
	virtual void Create_SC_w_Mat( SparseMatrix & K_in, SparseMatrix & B_in, SparseMatrix & SC_out, bool isThreaded, MKL_INT generate_symmetric_sc_1_generate_general_sc_0 ) = 0;
	virtual void Create_non_sym_SC_w_Mat( SparseMatrix & K_in, SparseMatrix & B1_in, SparseMatrix & B0_in, SparseMatrix & SC_out, bool isThreaded, MKL_INT generate_symmetric_sc_1_generate_general_sc_0 ) = 0;

	virtual void SolveCG(SparseMatrix & A_in, SEQ_VECTOR <double> & rhs_in, SEQ_VECTOR <double> & sol, SEQ_VECTOR <double> & initial_guess) = 0;
	virtual void SolveCG(SparseMatrix & A_in, SEQ_VECTOR <double> & rhs, SEQ_VECTOR <double> & sol) = 0;
	virtual void SolveCG(SparseMatrix & A_in, SEQ_VECTOR <double> & rhs_sol) = 0;
};



#endif /* SOLVER_SPECIFIC_SPARSESOLVER_H_ */
