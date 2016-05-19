#ifndef SPARSE_SOLVER_CUDA_H_
#define SPARSE_SOLVER_CUDA_H_

#include <cusolverSp.h>
#include <cusolverSp_LOWLEVEL_PREVIEW.h>

#include "../sparsesolver.h"

using std::string;
using std::endl;
using std::left;
using std::fixed;

#include <dmumps_c.h>


#pragma once

namespace espreso {

class SparseSolverCUDA: public SparseSolver
{

public:
	// Constructor
	SparseSolverCUDA();

	//Destructor
	~SparseSolverCUDA();

	DMUMPS_STRUC_C id;

	bool 		initialized;
	bool 		keep_factors;
	bool 		import_with_copy;
	int  		MPIrank;
	bool 		USE_FLOAT;
	int  		reorder; // 0 - no reordering, 1 - symrcm, 2 - symamd

	// Matrix properties
	MKL_INT 	rows;
	MKL_INT 	cols;
	MKL_INT 	nnz;

	MKL_INT		I_row_indices_size;
	MKL_INT		J_col_indices_size;
	MKL_INT		V_values_size;

	// MKL_INT		* I_row_indices;
	// MKL_INT		* J_col_indices;
	// double		* V_values;

	// MKL_INT		* CSR_I_row_indices;
	// MKL_INT		* CSR_J_col_indices;
	// double		* CSR_V_values;
	// float		* CSR_V_values_fl;

	int*		permutation;

	MKL_INT		CSR_I_row_indices_size;
	MKL_INT		CSR_J_col_indices_size;
	MKL_INT		CSR_V_values_size;
	MKL_INT		CSR_V_values_fl_size;
	MKL_INT		permutation_size;

	SEQ_VECTOR <double>	rhs_sol_reordered;	// only reordered
	SEQ_VECTOR <float>	rhs_sol_fl; 	// original or reordered

	// cuSolver variables
	cudaStream_t		cuStream;
	cusolverSpHandle_t	soHandle;
	csrcholInfo_t		soInfo;
	cusparseMatDescr_t	matDescr;

	// CUDA device variables
	int*		D_CSR_I_row_indices;
	int*		D_CSR_J_col_indices;
	double*		D_CSR_V_values;
	double*		D_rhs_sol;
	float*		D_CSR_V_values_fl;
	float*		D_rhs_sol_fl;
	void*		D_buffer;

	size_t 		internalDataInBytes;
	size_t 		workspaceInBytes;
	bool		keep_buffer;

	// *** Pardiso Solver Variables
	MKL_INT 	mtype;		/* Real symmetric matrix */
	MKL_INT 	iparm[65]; // typ matice
	MKL_INT 	maxfct, mnum, phase, error, msglvl;
	// ***

	MKL_INT 	m_nRhs;
	MKL_INT 	m_factorized;
	MKL_INT 	m_Kplus_size;
	// END - MKL DSS Solver Variables

	// Matrices
	//SparseMatrix m_A;

	// for in-place solve
	// SEQ_VECTOR <double> tmp_sol;

	//Members
	void ReorderMatrix(SparseMatrix & A);

	void ImportMatrix(SparseMatrix & A);
	void ImportMatrix_fl(SparseMatrix & A);
	void ImportMatrix_wo_Copy(SparseMatrix & A);

	int Factorization(const std::string &str);
	void Clear();
	void SetThreaded();


	void Solve( SEQ_VECTOR <double> & rhs, SEQ_VECTOR <double> & sol, MKL_INT rhs_start_index, MKL_INT sol_start_index);
	void Solve( SEQ_VECTOR <double> & rhs, SEQ_VECTOR <double> & sol, MKL_INT n_rhs );
	void Solve( SEQ_VECTOR <double> & rhs_sol);

	void SolveMat_Dense( SparseMatrix & A_in_out );
	void SolveMat_Dense( SparseMatrix & A_in, SparseMatrix & B_out );

	void SolveMatF( SparseMatrix & A_in, SparseMatrix & B_out, bool isThreaded );

	void SolveMat_Sparse( SparseMatrix & A );
	void SolveMat_Sparse( SparseMatrix & A_in, SparseMatrix & B_out );
	void SolveMat_Sparse( SparseMatrix & A_in, SparseMatrix & B_out, char T_for_input_matrix_is_transposed_N_input_matrix_is_NOT_transposed );

	void Create_SC( SparseMatrix & B_out, MKL_INT sc_size, bool isThreaded );
	void Create_SC_w_Mat( SparseMatrix & K_in, SparseMatrix & B_in, SparseMatrix & SC_out, bool isThreaded, MKL_INT generate_symmetric_sc_1_generate_general_sc_0 );
	void Create_non_sym_SC_w_Mat( SparseMatrix & K_in, SparseMatrix & B1_in, SparseMatrix & B0_in, SparseMatrix & SC_out, bool isThreaded, MKL_INT generate_symmetric_sc_1_generate_general_sc_0 );

	void SolveCG(SparseMatrix & A_in, SEQ_VECTOR <double> & rhs_in, SEQ_VECTOR <double> & sol, SEQ_VECTOR <double> & initial_guess);
	void SolveCG(SparseMatrix & A_in, SEQ_VECTOR <double> & rhs, SEQ_VECTOR <double> & sol);
	void SolveCG(SparseMatrix & A_in, SEQ_VECTOR <double> & rhs_sol);
};

}

#endif