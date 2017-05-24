#ifndef DENSE_SOLVER_CUDA_H_
#define DENSE_SOLVER_CUDA_H_

#include <cusolverDn.h>
#include <cublas_v2.h>

#include "../densesolver.h"

// using std::string;
// using std::endl;
// using std::left;
// using std::fixed;

// #include <dmumps_c.h>


#pragma once

namespace espreso {

class DenseSolverCUDA: public DenseSolver
{

public:
	// Constructor
	DenseSolverCUDA();

	//Destructor
	~DenseSolverCUDA();

	// DMUMPS_STRUC_C id;

	// bool 		initialized;
	// bool 		keep_factors;
	bool 		import_with_copy;
	// int  		MPIrank;
	bool 		USE_FLOAT;

	// Matrix properties
	MKL_INT 	m_rows;
	MKL_INT 	m_cols;
	MKL_INT 	m_nnz;
	MKL_INT		m_lda;
	MKL_INT		m_ldb;

	// Dense data
	MKL_INT		m_dense_values_size;
	MKL_INT		m_dense_values_fl_size;
	float *  	m_dense_values_fl;
	
	SEQ_VECTOR <float>		tmp_sol_fl;

	// cuSolver variables
	cudaStream_t		cuStream;
	cusolverDnHandle_t 	soDnHandle;

	// CUDA device variables
	int * 		D_devInfo;
	double * 	D_dense_values;
	float * 	D_dense_values_fl;
	double *	D_B_dense_values;
	float *		D_B_dense_values_fl;
	
	// double*		D_rhs_sol;
	// float*		D_rhs_sol_fl;

	// bool			keep_buffer;

	// *** Pardiso Solver Variables
	// MKL_INT 	mtype;		/* Real symmetric matrix */
	MKL_INT 	iparm[65]; // typ matice
	// MKL_INT 	maxfct, mnum, phase, error;
	// ***

	MKL_INT 	m_nRhs;
	// MKL_INT 	m_factorized;

	// Matrices
	// //SparseMatrix m_A;

	// // SEQ_VECTOR <double> tmp_sol;

	//Members
	void ImportMatrix(SparseMatrix & A);
	void ImportMatrix_fl(SparseMatrix & A);
	void ImportMatrix_wo_Copy(SparseMatrix & A);

	void ImportMatrix_wo_Copy_fl(SparseMatrix & A) {};

	int Factorization(const std::string &str);
	void Clear();
	void SetThreaded();


	void Solve( SEQ_VECTOR <double> & rhs, SEQ_VECTOR <double> & sol, MKL_INT rhs_start_index, MKL_INT sol_start_index);
	void Solve( SEQ_VECTOR <double> & rhs, SEQ_VECTOR <double> & sol, MKL_INT n_rhs );
	void Solve( SEQ_VECTOR <double> & rhs_sol);

	void SolveMat_Sparse( SparseMatrix & A );
	void SolveMat_Sparse( SparseMatrix & A_in, SparseMatrix & B_out );
	void SolveMat_Sparse( SparseMatrix & A_in, SparseMatrix & B_out, char T_for_input_matrix_is_transposed_N_input_matrix_is_NOT_transposed );

	void SolveMat_Dense( SparseMatrix & A_in_out );
	void SolveMat_Dense( SparseMatrix & A_in, SparseMatrix & B_out );

	void SolveCG(SparseMatrix & A_in, SEQ_VECTOR <double> & rhs_in, SEQ_VECTOR <double> & sol, SEQ_VECTOR <double> & initial_guess);
	void SolveCG(SparseMatrix & A_in, std::vector <double> & rhs, std::vector <double> & sol);
	void SolveCG(SparseMatrix & A_in, std::vector <double> & rhs_sol);
};

}

#endif
