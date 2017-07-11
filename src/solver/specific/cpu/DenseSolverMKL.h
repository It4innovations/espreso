#ifndef DENSE_SOLVER_MKL_H_
#define DENSE_SOLVER_MKL_H_

//#include "mkl.h"

#include "mkl_lapack.h"

#include "../densesolver.h"

// using std::string;
// using std::endl;
// using std::left;
// using std::fixed;

// #include <dmumps_c.h>


#pragma once

namespace espreso {

class DenseSolverMKL: public DenseSolver
{

public:
	// Constructor
	DenseSolverMKL();

	//Destructor
	~DenseSolverMKL();

	// DMUMPS_STRUC_C id;

	eslocal m_Kplus_size;
	bool keep_factors;
	eslocal msglvl;
	eslocal MPIrank;



	// bool 		initialized;
	// bool 		keep_factors;
	bool 		import_with_copy;
	// int  		MPIrank;
	bool 		USE_FLOAT;

	//eslocal           mtype;
	MatrixType        mtype;
	//MKL_INT mtype;		/* Real symmetric matrix */


	// Matrix properties
	MKL_INT 	m_rows;
	MKL_INT 	m_cols;
	MKL_INT 	m_nnz;
	MKL_INT		m_lda;

	// Dense data
	MKL_INT					m_dense_values_size;
	MKL_INT					m_dense_values_fl_size;
	double * 				m_dense_values;
	float *  				m_dense_values_fl;
	SEQ_VECTOR <eslocal>    m_ipiv;

	SEQ_VECTOR <float>		tmp_sol_fl;

	// *** Pardiso Solver Variables
	// MKL_INT 	mtype;		/* Real symmetric matrix */
	MKL_INT 	iparm[65]; // typ matice
	// MKL_INT 	maxfct, mnum, phase, error, msglvl;
	// ***

	MKL_INT 	m_nRhs;
	// MKL_INT 	m_factorized;

	//Members
	void ImportMatrix(SparseMatrix & A);
	void ImportMatrix_fl(SparseMatrix & A);
	void ImportMatrix_wo_Copy(SparseMatrix & A);

	void ImportMatrix_wo_Copy_fl(SparseMatrix & A) {};

	int  Factorization(const std::string &str);
	void Clear();
	void SetThreaded();


	void Solve( SEQ_VECTOR <double> & rhs, SEQ_VECTOR <double> & sol, MKL_INT rhs_start_index, MKL_INT sol_start_index);
	void Solve( SEQ_VECTOR <double> & rhs, SEQ_VECTOR <double> & sol, MKL_INT n_rhs );
	void Solve( SEQ_VECTOR <double> & rhs_sol, MKL_INT n_rhs);
	void Solve( SEQ_VECTOR <double> & rhs_sol);

	void SolveMat_Sparse( SparseMatrix & A );
	void SolveMat_Sparse( SparseMatrix & A_in, SparseMatrix & B_out );
	void SolveMat_Sparse( SparseMatrix & A_in, SparseMatrix & B_out, char T_for_input_matrix_is_transposed_N_input_matrix_is_NOT_transposed );

	void SolveMat_Dense( SparseMatrix & A_in_out );
	void SolveMat_Dense( SparseMatrix & A_in, SparseMatrix & B_out );

	void SolveCG(SparseMatrix & A_in, SEQ_VECTOR <double> & rhs_in, SEQ_VECTOR <double> & sol, SEQ_VECTOR <double> & initial_guess){};
	void SolveCG(SparseMatrix & A_in, std::vector <double> & rhs, std::vector <double> & sol){};
	void SolveCG(SparseMatrix & A_in, std::vector <double> & rhs_sol){};


};

}

#endif
