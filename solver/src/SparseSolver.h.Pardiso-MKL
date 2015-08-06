#ifdef WIN32	 
	#include "stdafx.h"
#endif

 

#include "utils.h"

//#include "mkl.h"
//#include <string>
//#include <sstream>
//#include <iostream>
//#include <vector>
//#include <fstream>
//#include <algorithm>
//
//#include <math.h>
//#include <stack>
//#include <ctime>

#include "SparseMatrix.h"

using std::string; 
using std::endl; 
using std::left;
using std::fixed;



#pragma once

class SparseSolver
{

public: 
	// Constructor 
	SparseSolver(); 
	
	//Destructor
	~SparseSolver();

	// Matrix properties 
	int rows; 
	int cols; 
	int nnz; 

	int		* CSR_I_row_indices;
	int		* CSR_J_col_indices;
	double	* CSR_V_values; 

	int		CSR_I_row_indices_size;
	int		CSR_J_col_indices_size;
	int		CSR_V_values_size; 

	// *** Pardiso Solver Variables 

	MKL_INT mtype;		/* Real symmetric matrix */

	/* Internal solver memory pointer pt, */
	/* 32-bit: int pt[64]; 64-bit: long int pt[64] */
	/* or void *pt[64] should be OK on both architectures */
	void *pt[64];
	//MKL_INT * pt; 


	/* Pardiso control parameters. */
	MKL_INT iparm[64];
	//MKL_INT * iparm;
	
	
	MKL_INT maxfct, mnum, phase, error, msglvl;

	// *** 

	int m_nRhs;
	int m_factorized; 
	int m_Kplus_size;
	// END - MKL DSS Solver Variables 

	// Matrices
	//SparseMatrix m_A; 

	// for in-place solve 
	SEQ_VECTOR <double> tmp_sol; 

	//Members
	void ImportMatrix(SparseMatrix & A);
	void Factorization(); 
	void Clear();
	
	void Solve( SEQ_VECTOR <double> & rhs, SEQ_VECTOR <double> & sol, int rhs_start_index, int sol_start_index);
	void Solve( SEQ_VECTOR <double> & rhs, SEQ_VECTOR <double> & sol, int n_rhs );
	void Solve( SEQ_VECTOR <double> & rhs_sol); 

	void SolveMat_Dense( SparseMatrix & A_in_out );
	void SolveMat_Dense( SparseMatrix & A_in, SparseMatrix & B_out ); 

	void SolveMatF( SparseMatrix & A_in, SparseMatrix & B_out ); 

	void SolveMat_Sparse( SparseMatrix & A ); 
	void SolveMat_Sparse( SparseMatrix & A_in, SparseMatrix & B_out ); 
	void SolveMat_Sparse( SparseMatrix & A_in, SparseMatrix & B_out, char T_for_input_matrix_is_transposed_N_input_matrix_is_NOT_transposed );

};