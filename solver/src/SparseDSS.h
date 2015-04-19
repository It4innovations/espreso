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

// *******************************************************************
// **** SparseDSS - Sparse Direct Solver Class ************************

class SparseDSS {

public: 
	// Constructor 
	SparseDSS(); 

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

	// MKL DSS Solver Variables 
	_MKL_DSS_HANDLE_t m_handle;
	_INTEGER_t m_error;

	MKL_INT m_opt;
	MKL_INT m_sym;
	MKL_INT m_type;
	MKL_INT	m_order;

	int m_nRhs;
	int m_factorized; 
	int m_Kplus_size;
	// END - MKL DSS Solver Variables 

	// Matrices
	SparseMatrix m_A; 


	//Members
	void LoadMatrix(string filename); 
	void LoadMatrixBin(string filename); 
	void ImportMatrix(SparseMatrix & A);
	void Factorization(); 
	void Solve( SEQ_VECTOR <double> & rhs, SEQ_VECTOR <double> & sol, int rhs_start_index, int sol_start_index);
	void Solve( SEQ_VECTOR <double> & rhs, SEQ_VECTOR <double> & sol, int n_rhs );

	void SolveMat_Dense( SparseMatrix & A_in_out );
	void SolveMat_Dense( SparseMatrix & A_in, SparseMatrix & B_out ); 

	void SolveMat_Sparse( SparseMatrix & A ); 
	void SolveMat_Sparse( SparseMatrix & A_in, SparseMatrix & B_out ); 
	void SolveMat_Sparse( SparseMatrix & A_in, SparseMatrix & B_out, char T_for_input_matrix_is_transposed_N_input_matrix_is_NOT_transposed );

};