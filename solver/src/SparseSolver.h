//#ifndef SPARSE_SOLVER_H_
//#define SPARSE_SOLVER_H_

#ifdef WIN32	 
	#include "stdafx.h"
#endif

 
#include "utils.h"
#include "SparseMatrix.h"

using std::string; 
using std::endl; 
using std::left;
using std::fixed;

#include <dmumps_c.h>


#pragma once

class SparseSolver
{

public: 
	// Constructor 
	SparseSolver(); 
	
	//Destructor
	~SparseSolver();

	DMUMPS_STRUC_C id;



	// Matrix properties 

	bool initialized;
	bool keep_factors;
	bool import_with_copy;
	int  MPIrank;

	MKL_INT rows;
	MKL_INT cols;
	MKL_INT nnz;



	MKL_INT		* I_row_indices;
	MKL_INT		* J_col_indices;
	double		* V_values;

	MKL_INT		I_row_indices_size;
	MKL_INT		J_col_indices_size;
	MKL_INT		V_values_size;


	MKL_INT		* CSR_I_row_indices;
	MKL_INT		* CSR_J_col_indices;
	double		* CSR_V_values;

	MKL_INT		CSR_I_row_indices_size;
	MKL_INT		CSR_J_col_indices_size;
	MKL_INT		CSR_V_values_size;

	// *** Pardiso Solver Variables 

	MKL_INT mtype;		/* Real symmetric matrix */

	/* Internal solver memory pointer pt, */
	/* 32-bit: int pt[64]; 64-bit: long int pt[64] */
	/* or void *pt[64] should be OK on both architectures */
	void *pt[64];
	//int * pt;


	/* Pardiso control parameters. */
	MKL_INT iparm[65];
	double  dparm[65];
	//int * iparm;
	
	MKL_INT maxfct, mnum, phase, error, msglvl;
	// *** 

	MKL_INT m_nRhs;
	MKL_INT m_factorized;
	MKL_INT m_Kplus_size;
	// END - MKL DSS Solver Variables 

	// Matrices
	//SparseMatrix m_A; 

	// for in-place solve 
	SEQ_VECTOR <double> tmp_sol; 

	//Members
	void ImportMatrix(SparseMatrix & A);
	void ImportMatrix_wo_Copy(SparseMatrix & A);

	void Factorization(); 
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

};
//#endif //SPARSE_SOLVER_H_
