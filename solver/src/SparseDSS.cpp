#include "SparseDSS.h"


SparseDSS::SparseDSS(){

	m_opt		 = MKL_DSS_DEFAULTS;
	m_sym		 = MKL_DSS_SYMMETRIC;
	m_type		 = MKL_DSS_POSITIVE_DEFINITE;
	m_order		 = MKL_DSS_AUTO_ORDER; 

	m_nRhs		 = 1; 
	m_factorized = 0; 

}

void SparseDSS::ImportMatrix(SparseMatrix & A) {

	rows	= A.rows;
	cols	= A.cols;
	nnz		= A.nnz; 
	m_Kplus_size = A.rows; 

	CSR_I_row_indices_size = A.CSR_I_row_indices.size();
	CSR_J_col_indices_size = A.CSR_J_col_indices.size(); 
	CSR_V_values_size	   = A.CSR_V_values.size(); 

	CSR_I_row_indices = (int*)malloc(CSR_I_row_indices_size * sizeof(int)); 
	CSR_J_col_indices = (int*)malloc(CSR_J_col_indices_size * sizeof(int));
	CSR_V_values =		(double*)malloc(CSR_V_values_size   * sizeof(double));

	copy(A.CSR_I_row_indices.begin(), A.CSR_I_row_indices.end(), CSR_I_row_indices); 
	copy(A.CSR_J_col_indices.begin(), A.CSR_J_col_indices.end(), CSR_J_col_indices); 
	copy(A.CSR_V_values.begin(),      A.CSR_V_values.end(),      CSR_V_values); 

}

void SparseDSS::LoadMatrix(string filename) {
	m_A.LoadMatrixBin(filename, 'S'); 
	ImportMatrix(m_A); 
	m_Kplus_size = m_A.rows;
	m_A.Clear(); 
}

void SparseDSS::LoadMatrixBin(string filename) {
	m_A.LoadMatrixBin(filename, 'S'); 
	ImportMatrix(m_A); 
	m_Kplus_size = m_A.rows;
	m_A.Clear(); 
}

void SparseDSS::Factorization() {

	/* --------------------- */
	/* Initialize the solver */
	/* --------------------- */
	m_error = dss_create (m_handle, m_opt);
	//if (error != MKL_DSS_SUCCESS)

	/* ------------------------------------------- */
	/* Define the non-zero structure of the matrix */
	/* ------------------------------------------- */
	//error = dss_define_structure (handle, sym, rowIndex,                nRows,   nCols,   columns,                 nNonZeros);
	//m_error =   dss_define_structure (m_handle, m_sym, &m_A.CSR_I_row_indices[0], m_A.rows,  m_A.cols,  &m_A.CSR_J_col_indices[0], m_A.nnz);
	m_error =   dss_define_structure (m_handle, m_sym, CSR_I_row_indices, rows,  cols,  CSR_J_col_indices, nnz);
	//if (error != MKL_DSS_SUCCESS)

	/* ------------------ */
	/* Reorder the matrix */
	/* ------------------ */
	//error = dss_reorder (handle, opt, 0);
	m_error = dss_reorder (m_handle, m_order, 0); 

	//if (error != MKL_DSS_SUCCESS)

	/* ------------------ */
	/* Factor the matrix */
	/* ------------------ */
	//error = dss_factor_real (handle, type, values);
	//m_error = dss_factor_real (m_handle, m_type, &m_A.CSR_V_values[0]);
	m_error = dss_factor_real (m_handle, m_type, CSR_V_values);
	//if (error != MKL_DSS_SUCCESS)
	m_factorized = 1; 
}

void SparseDSS::Solve( SEQ_VECTOR <double> & rhs, SEQ_VECTOR <double> & sol, int rhs_start_index, int sol_start_index) {

	int nRhs = 1; 

	// if Matrix is not factorized yet, run factorization
	//	if (m_factorized == 0)
	//		Factorization(); 

	// if sol vector size is small, resize it 
	//if (sol.size() < A.rows)
	//	sol.resize(A.rows);

	/* ------------------------ */
	/* Get the solution vector */
	/* ------------------------ */
	//m_error = dss_solve_real (m_handle, m_opt, &rhs[rhs_start_index], m_nRhs, &sol[sol_start_index]);
	m_error = dss_solve_real (m_handle, m_opt, &rhs[rhs_start_index], nRhs, &sol[sol_start_index]);
	//if (error != MKL_DSS_SUCCESS)

	/* -------------------------- */
	/* Deallocate solver storage */
	/* -------------------------- */
	//m_error = dss_delete (m_handle, m_opt);
	//if (error != MKL_DSS_SUCCESS) 
} 

void SparseDSS::Solve( SEQ_VECTOR <double> & rhs, SEQ_VECTOR <double> & sol, int n_rhs) {
	
	//if (m_factorized == 0)
	//	Factorization(); 

	m_error = dss_solve_real (m_handle, m_opt, &rhs[0], n_rhs, &sol[0]);
}

void SparseDSS::SolveMat_Dense( SparseMatrix & A ) {
	SolveMat_Dense(A, A);
}

void SparseDSS::SolveMat_Dense( SparseMatrix & A_in, SparseMatrix & B_out ) {

	//if (m_factorized == 0)
	//	Factorization(); 

	SEQ_VECTOR<double> rhs; 
	SEQ_VECTOR<double> sol;

	int job[8];
	job[0] = 1; // if job(1)=1, the rectangular matrix A is restored from the CSR format.
	job[1] = 1; // if job(2)=1, one-based indexing for the rectangular matrix A is used.
	job[2] = 1; // if job(3)=1, one-based indexing for the matrix in CSR format is used.
	job[3] = 2; // If job(4)=2, adns is a whole matrix A.

	job[4] = 0; // job(5)=nzmax - maximum number of the non-zero elements allowed if job(1)=0.
	job[5] = 0; // job(6) - job indicator for conversion to CSR format.
	// If job(6)=0, only array ia is generated for the output storage.
	// If job(6)>0, arrays acsr, ia, ja are generated for the output storage.
	job[6] = 0; //
	job[7] = 0; //

	MKL_INT m = A_in.rows;
	MKL_INT n = A_in.cols; 
	MKL_INT nRhs = A_in.cols;

	rhs.resize(m * n);
	sol.resize(m * n);
	//double *Adns = &rhs[0]; 
	MKL_INT lda  = m;

	//double *Acsr = &A_in.CSR_V_values[0];
	//MKL_INT *AJ  = &A_in.CSR_J_col_indices[0];
	//MKL_INT *AI  = &A_in.CSR_I_row_indices[0];

	MKL_INT info;


	// Convert input matrix (RHS) to dense format 

	//void mkl_ddnscsr (
	//	MKL_INT *job, 
	//	MKL_INT *m, MKL_INT *n, 
	//	double *Adns, MKL_INT *lda, 
	//	double *Acsr, MKL_INT *AJ, MKL_INT *AI, 
	//	MKL_INT *info);

	mkl_ddnscsr (
		job, 
		&m, &n, 
		&rhs[0], &lda, 
		&A_in.CSR_V_values[0], &A_in.CSR_J_col_indices[0], &A_in.CSR_I_row_indices[0], 
		&info);

	// Solve with multiple right hand sides 
	m_error = dss_solve_real (m_handle, m_opt, &rhs[0], nRhs, &sol[0]);
	rhs.clear();

	// Convert solution matrix (SOL) to sparse format - find nnz step
	job[0] = 0; // If job(1)=0, the rectangular matrix A is converted to the CSR format;
	job[1] = 1; // if job(2)=1, one-based indexing for the rectangular matrix A is used.
	job[2] = 1; // if job(3)=1, one-based indexing for the matrix in CSR format is used.
	job[3] = 2; // If job(4)=2, adns is a whole matrix A.

	job[4] = 1; // job(5)=nzmax - maximum number of the non-zero elements allowed if job(1)=0.
	job[5] = 0; // job(6) - job indicator for conversion to CSR format.
	// If job(6)=0, only array ia is generated for the output storage.
	// If job(6)>0, arrays acsr, ia, ja are generated for the output storage.
	job[6] = 0; //
	job[7] = 0; //

	//Adns = &sol[0]; 

	B_out.CSR_I_row_indices.resize(m + 1); 
	B_out.CSR_J_col_indices.resize(1);
	B_out.CSR_V_values.resize(1);

	//Acsr = &B_out.CSR_V_values[0];
	//AJ   = &B_out.CSR_J_col_indices[0];
	//AI   = &B_out.CSR_I_row_indices[0];

	mkl_ddnscsr (
		job, 
		&m, &n, 
		&sol[0], &lda, 
		&B_out.CSR_V_values[0], &B_out.CSR_J_col_indices[0], &B_out.CSR_I_row_indices[0], 
		&info);

	// Convert solution matrix (SOL) to sparse format - convert step 
	int nnzmax = B_out.CSR_I_row_indices[m]-1; 

	B_out.CSR_J_col_indices.resize(nnzmax);
	B_out.CSR_V_values.resize(nnzmax);

	//Acsr = &B_out.CSR_V_values[0];
	//AJ   = &B_out.CSR_V_values;
	//AI   = &B_out.CSR_I_row_indices[0];

	job[4] = nnzmax; // job(5) = nzmax - maximum number of the non-zero elements allowed if job(1)=0.
	job[5] = 1; // job(6) - job indicator for conversion to CSR format.
	// If job(6)=0, only array ia is generated for the output storage.
	// If job(6)>0, arrays acsr, ia, ja are generated for the output storage.

	mkl_ddnscsr (
		job, 
		&m, &n, 
		&sol[0], &lda, 
		&B_out.CSR_V_values[0], &B_out.CSR_J_col_indices[0], &B_out.CSR_I_row_indices[0], 
		&info);


	// Setup parameters for output matrix 
	B_out.cols	= A_in.cols;
	B_out.rows	= A_in.rows;
	B_out.nnz	= B_out.CSR_V_values.size();
	B_out.type	= 'G';

	sol.clear(); 

}

void SparseDSS::SolveMat_Sparse( SparseMatrix & A) {
	SolveMat_Sparse(A, A); 
};

void SparseDSS::SolveMat_Sparse( SparseMatrix & A_in, SparseMatrix & B_out) {
	SolveMat_Sparse(A_in, B_out, 'T'); 
};

void SparseDSS::SolveMat_Sparse( SparseMatrix & A_in, SparseMatrix & B_out, char T_for_input_matrix_is_transposed_N_input_matrix_is_NOT_transposed ) {

	char trans = T_for_input_matrix_is_transposed_N_input_matrix_is_NOT_transposed;

	SparseMatrix tmpM;
	if (trans == 'T')
		A_in.MatTranspose(tmpM);
	else
		tmpM = A_in; 

	//if (m_factorized == 0)
	//	Factorization(); 

	SEQ_VECTOR<double> rhs; 
	SEQ_VECTOR<double> sol;

	rhs.resize(tmpM.cols);
	sol.resize(tmpM.cols);

	// main loop over rows 
	int col = 0; 
	int n_nnz = 0; 
	for (int row = 1; row < tmpM.CSR_I_row_indices.size(); row++) {
		int row_size = tmpM.CSR_I_row_indices[row] - tmpM.CSR_I_row_indices[row-1];
		if (row_size > 0) {
			for (int c = 0; c < row_size; c++) { // loop over selected row 
				rhs[ tmpM.CSR_J_col_indices[col] - 1] = tmpM.CSR_V_values [col];
				col++;
			}
			int nRhs_l = 1;
			m_error = dss_solve_real (m_handle, m_opt, &rhs[0], nRhs_l, &sol[0]);

			for (int s = 0; s < sol.size(); s++){
				if (sol[s] != 0.0) {
					tmpM.I_row_indices.push_back(row);
					tmpM.J_col_indices.push_back(s+1);
					tmpM.V_values.push_back(sol[s]);
					n_nnz++; 
				}
			}

			//Reset RHS and SOL
			fill(rhs.begin(), rhs.end(), 0); // reset entire vector to 0
			//fill(sol.begin(), sol.end(), 0); // reset entire vector to 0
		}
	}

	rhs.clear();
	sol.clear();

	tmpM.nnz = n_nnz; 
	tmpM.ConvertToCSR(1); 
	tmpM.MatTranspose(B_out); 

	tmpM.Clear(); 

}

// **** END - SparseDS - Sparse Direct Solver Class*******************
// *******************************************************************