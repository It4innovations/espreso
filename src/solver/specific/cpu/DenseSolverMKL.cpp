
#include "DenseSolverMKL.h"

using namespace espreso;

DenseSolverMKL::DenseSolverMKL(){

	// keep_factors=true;
	// initialized = false;
	USE_FLOAT = false;
	import_with_copy = false;
	mtype = espreso::SparseMatrix::MatrixType::REAL_UNSYMMETRIC;

	m_dense_values_size = 0;
	m_dense_values_fl_size = 0;

	/* Numbers of processors, value of OMP_NUM_THREADS */
	int num_procs;
	char * var = getenv("PAR_NUM_THREADS");
   	if(var != NULL)
   		sscanf( var, "%d", &num_procs );
	else {
   		printf("Set environment PAR_NUM_THREADS to 1");
    	exit(1);
	}

    iparm[2]  = num_procs;

	m_nRhs		 = 1;
// 	m_factorized = 0;
}

DenseSolverMKL::~DenseSolverMKL() {

	this->Clear();
}

void DenseSolverMKL::Clear() {

	if (import_with_copy) {
		if(m_dense_values_size > 0)	delete [] m_dense_values;
		if(m_dense_values_fl_size > 0)	delete [] m_dense_values_fl;
	}

	m_dense_values_size = 0;
	m_dense_values_fl_size = 0;

	m_dense_values = NULL;
	m_dense_values_fl = NULL;
}

void DenseSolverMKL::ImportMatrix(espreso::SparseMatrix & A) {

	USE_FLOAT 	= false;
	
	m_rows		= A.rows;
	m_cols		= A.cols;
	m_nnz		= A.nnz;
	m_lda		= A.rows;

	m_dense_values_size = A.dense_values.size();
	m_dense_values_fl_size = 0;

	m_dense_values 	= new double[m_dense_values_size];

	copy(A.dense_values.begin(), A.dense_values.end(), m_dense_values);

	import_with_copy = true;

	mtype = A.mtype;
}

void DenseSolverMKL::ImportMatrix_fl(espreso::SparseMatrix & A) {

	USE_FLOAT 	= true;

	m_rows		= A.rows;
	m_cols		= A.cols;
	m_nnz		= A.nnz;
	m_lda		= A.rows;

	m_dense_values_size = 0;
	m_dense_values_fl_size = A.dense_values.size();

	m_dense_values_fl = new float[m_dense_values_fl_size];

	for (eslocal i = 0; i < m_dense_values_fl_size; i++)
		m_dense_values_fl[i] = (float) A.dense_values[i];

	import_with_copy = true;

	mtype = A.mtype;
}


void DenseSolverMKL::ImportMatrix_wo_Copy(espreso::SparseMatrix & A) {

	USE_FLOAT = false;

	m_rows		= A.rows;
	m_cols		= A.cols;
	m_nnz		= A.nnz;
	m_lda		= A.rows;

	m_dense_values_size = A.dense_values.size();
	m_dense_values_fl_size = 0;

	m_dense_values = &A.dense_values[0];
	
	import_with_copy = false;

	mtype = A.mtype;

}

void DenseSolverMKL::SetThreaded() {

	/* Numbers of processors, value of OMP_NUM_THREADS */
	int num_procs = Esutils::getEnv<int>("SOLVER_NUM_THREADS");

    iparm[2]  = num_procs;
}

int DenseSolverMKL::Factorization(const std::string &str) {

	eslocal info;
	char U = 'U';

	m_ipiv.resize(m_cols);

	switch (mtype) {
	case espreso::SparseMatrix::MatrixType::REAL_SYMMETRIC_POSITIVE_DEFINITE:

		if (USE_FLOAT) {
			ssptrf( &U, &m_cols, &m_dense_values_fl[0], &m_ipiv[0] , &info );
		} else {
			dsptrf( &U, &m_cols, &m_dense_values[0], &m_ipiv[0] , &info );
		}

		break;
	case espreso::SparseMatrix::MatrixType::REAL_SYMMETRIC_INDEFINITE:

		if (USE_FLOAT) {
			ssptrf( &U, &m_cols, &m_dense_values_fl[0], &m_ipiv[0] , &info );
		} else {
			dsptrf( &U, &m_cols, &m_dense_values[0], &m_ipiv[0] , &info );
		}

		break;
	case espreso::SparseMatrix::MatrixType::REAL_UNSYMMETRIC:

		if (USE_FLOAT) {
			//sgetrf( m, n, a, lda, ipiv, info )
			sgetrf( &m_cols, &m_cols, &m_dense_values_fl[0], &m_cols, &m_ipiv[0] , &info );
		} else {
			dgetrf( &m_cols, &m_cols, &m_dense_values[0], &m_cols, &m_ipiv[0] , &info );
		}

		break;
	}


	return info;
}

void DenseSolverMKL::Solve( SEQ_VECTOR <double> & rhs_sol) {

	char U = 'U';
	eslocal info = 0;
	m_nRhs = 1;



	switch (mtype) {
	case espreso::SparseMatrix::MatrixType::REAL_SYMMETRIC_POSITIVE_DEFINITE:
		if (USE_FLOAT) {
			tmp_sol_fl.resize(rhs_sol.size());

			for (size_t i = 0; i < rhs_sol.size(); i++)
				tmp_sol_fl[i] = (float)rhs_sol[i];

			ssptrs( &U, &m_rows, &m_nRhs, &m_dense_values_fl[0], &m_ipiv[0], &tmp_sol_fl[0], &m_rows, &info );

			for (size_t i = 0; i < rhs_sol.size(); i++)
				rhs_sol[i] = (double)tmp_sol_fl[i];
		} else {
			dsptrs( &U, &m_rows, &m_nRhs, &m_dense_values[0], &m_ipiv[0], &rhs_sol[0], &m_rows, &info );
		}

		break;
	case espreso::SparseMatrix::MatrixType::REAL_SYMMETRIC_INDEFINITE:

		if (USE_FLOAT) {
			tmp_sol_fl.resize(rhs_sol.size());

			for (size_t i = 0; i < rhs_sol.size(); i++)
				tmp_sol_fl[i] = (float)rhs_sol[i];

			ssptrs( &U, &m_rows, &m_nRhs, &m_dense_values_fl[0], &m_ipiv[0], &tmp_sol_fl[0], &m_rows, &info );

			for (size_t i = 0; i < rhs_sol.size(); i++)
				rhs_sol[i] = (double)tmp_sol_fl[i];
		} else {
			dsptrs( &U, &m_rows, &m_nRhs, &m_dense_values[0], &m_ipiv[0], &rhs_sol[0], &m_rows, &info );
		}

		break;
	case espreso::SparseMatrix::MatrixType::REAL_UNSYMMETRIC:

		if (USE_FLOAT) {
			tmp_sol_fl.resize(rhs_sol.size());

			for (size_t i = 0; i < rhs_sol.size(); i++)
				tmp_sol_fl[i] = (float)rhs_sol[i];

			ssptrs( &U, &m_rows, &m_nRhs, &m_dense_values_fl[0], &m_ipiv[0], &tmp_sol_fl[0], &m_rows, &info );
			//call sgetrs( trans, n,        nrhs,    a,                    lda,          ipiv,               b,     ldb, info )
			char trans = 'N';
			sgetrs(       &trans, &m_rows, &m_nRhs, &m_dense_values_fl[0], &m_rows, &m_ipiv[0], &tmp_sol_fl[0], &m_nRhs, &info);

			for (size_t i = 0; i < rhs_sol.size(); i++)
				rhs_sol[i] = (double)tmp_sol_fl[i];
		} else {
			//dsptrs( &U, &m_rows, &m_nRhs, &m_dense_values[0], &m_ipiv[0], &rhs_sol[0], &m_rows, &info );
			char trans = 'N';
			dgetrs(       &trans, &m_rows, &m_nRhs, &m_dense_values[0], &m_rows, &m_ipiv[0], &rhs_sol[0], &m_nRhs, &info);
		}

		break;
	}

}

void DenseSolverMKL::Solve( SEQ_VECTOR <double> & rhs, SEQ_VECTOR <double> & sol, MKL_INT n_rhs) {

	char U = 'U';
	eslocal info = 0;
	eslocal m_nRhs = n_rhs;
//
//	if (USE_FLOAT) {
//		sol.resize(rhs.size());
//		tmp_sol_fl.resize(rhs.size());
//
//		for (eslocal i = 0; i < rhs.size(); i++)
//			tmp_sol_fl[i] = (float)rhs[i];
//
//		ssptrs( &U, &m_rows, &n_rhs, &m_dense_values_fl[0], &m_ipiv[0], &tmp_sol_fl[0], &m_rows, &info );
//
//		for (size_t i = 0; i < rhs.size(); i++)
//			sol[i] = (double)tmp_sol_fl[i];
//	} else {
//		sol = rhs;
//		dsptrs( &U, &m_rows, &n_rhs, &m_dense_values[0], &m_ipiv[0], &sol[0], &m_rows, &info );
//	}

	switch (mtype) {
	case espreso::SparseMatrix::MatrixType::REAL_SYMMETRIC_POSITIVE_DEFINITE:
		if (USE_FLOAT) {
			sol.resize(rhs.size());
			tmp_sol_fl.resize(rhs.size());

			for (size_t i = 0; i < rhs.size(); i++)
				tmp_sol_fl[i] = (float)rhs[i];

			ssptrs( &U, &m_rows, &n_rhs, &m_dense_values_fl[0], &m_ipiv[0], &tmp_sol_fl[0], &m_rows, &info );

			for (size_t i = 0; i < rhs.size(); i++)
				sol[i] = (double)tmp_sol_fl[i];
		} else {
			sol = rhs;
			dsptrs( &U, &m_rows, &n_rhs, &m_dense_values[0], &m_ipiv[0], &sol[0], &m_rows, &info );
		}

		break;
	case espreso::SparseMatrix::MatrixType::REAL_SYMMETRIC_INDEFINITE:

		if (USE_FLOAT) {
			sol.resize(rhs.size());
			tmp_sol_fl.resize(rhs.size());

			for (size_t i = 0; i < rhs.size(); i++)
				tmp_sol_fl[i] = (float)rhs[i];

			ssptrs( &U, &m_rows, &n_rhs, &m_dense_values_fl[0], &m_ipiv[0], &tmp_sol_fl[0], &m_rows, &info );

			for (size_t i = 0; i < rhs.size(); i++)
				sol[i] = (double)tmp_sol_fl[i];
		} else {
			sol = rhs;
			dsptrs( &U, &m_rows, &n_rhs, &m_dense_values[0], &m_ipiv[0], &sol[0], &m_rows, &info );
		}

		break;
	case espreso::SparseMatrix::MatrixType::REAL_UNSYMMETRIC:

		char trans = 'N';

		if (USE_FLOAT) {
			sol.resize(rhs.size());
			tmp_sol_fl.resize(rhs.size());

			for (size_t i = 0; i < rhs.size(); i++)
				tmp_sol_fl[i] = (float)rhs[i];

			//ssptrs( &U, &m_rows, &n_rhs, &m_dense_values_fl[0], &m_ipiv[0], &tmp_sol_fl[0], &m_rows, &info );
			sgetrs(       &trans, &m_rows, &m_nRhs, &m_dense_values_fl[0], &m_rows, &m_ipiv[0], &tmp_sol_fl[0], &m_nRhs, &info);


			for (size_t i = 0; i < rhs.size(); i++)
				sol[i] = (double)tmp_sol_fl[i];
		} else {
			sol = rhs;
			//dsptrs( &U, &m_rows, &n_rhs, &m_dense_values[0], &m_ipiv[0], &sol[0], &m_rows, &info );
			//dgetrs(      trans,       n,    nrhs,                  a,     lda,       ipiv,       b,     ldb,  info )
			dgetrs(       &trans, &m_rows, &m_nRhs, &m_dense_values[0], &m_rows, &m_ipiv[0], &sol[0], &m_rows, &info);

		}

		break;
	}



}

void DenseSolverMKL::Solve( SEQ_VECTOR <double> & rhs, SEQ_VECTOR <double> & sol, MKL_INT rhs_start_index, MKL_INT sol_start_index) {

	printf("DenseSolverMKL::Solve( SEQ_VECTOR <double> & rhs, SEQ_VECTOR <double> & sol, MKL_INT rhs_start_index, MKL_INT sol_start_index) not implemented yet.\n");
	exit(1);

}



void DenseSolverMKL::SolveMat_Sparse( espreso::SparseMatrix & A) {
	SolveMat_Sparse(A, A);
};

void DenseSolverMKL::SolveMat_Sparse( espreso::SparseMatrix & A_in, espreso::SparseMatrix & B_out) {
	SolveMat_Sparse(A_in, B_out, 'T');
};

void DenseSolverMKL::SolveMat_Sparse( espreso::SparseMatrix & A_in, espreso::SparseMatrix & B_out, char T_for_input_matrix_is_transposed_N_input_matrix_is_NOT_transposed ) {

	char trans = T_for_input_matrix_is_transposed_N_input_matrix_is_NOT_transposed;

	espreso::SparseMatrix tmpM;
	if (trans == 'T')
		A_in.MatTranspose(tmpM);
	else
		tmpM = A_in;


	SEQ_VECTOR<double> rhs;
	SEQ_VECTOR<double> sol;

	rhs.reserve(tmpM.cols+8);
	rhs.reserve(tmpM.cols+8);

	rhs.resize(tmpM.cols);
	sol.resize(tmpM.cols);

	// main loop over rows
	MKL_INT col = 0;
	MKL_INT n_nnz = 0;
	for (size_t row = 1; row < tmpM.CSR_I_row_indices.size(); row++) {
		MKL_INT row_size = tmpM.CSR_I_row_indices[row] - tmpM.CSR_I_row_indices[row-1];
		if (row_size > 0) {
			for (MKL_INT c = 0; c < row_size; c++) { // loop over selected row
				rhs[ tmpM.CSR_J_col_indices[col] - 1] = tmpM.CSR_V_values [col];
				col++;
			}
			MKL_INT nRhs_l = 1;
			//m_error = dss_solve_real (m_handle, m_opt, &rhs[0], nRhs_l, &sol[0]);
			Solve(rhs, sol, nRhs_l);

			for (size_t s = 0; s < sol.size(); s++){
				if (sol[s] != 0.0) {
					tmpM.I_row_indices.push_back(row);
					tmpM.J_col_indices.push_back(s+1);
					tmpM.V_values.push_back(sol[s]);
					n_nnz++;
				}
			}

			//Reset InitialCondition and SOL
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

