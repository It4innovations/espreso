
#include "DenseSolverMKL.h"

using namespace espreso;

DenseSolverMKL::DenseSolverMKL(){

	// keep_factors=true;
	// initialized = false;
	USE_FLOAT = false;
	import_with_copy = false;

	m_dense_values_size = 0;
	m_dense_values_fl_size = 0;

	/* Numbers of processors, value of OMP_NUM_THREADS */
	int num_procs;
	char * var = getenv("OMP_NUM_THREADS");
   	if(var != NULL)
   		sscanf( var, "%d", &num_procs );
	else {
   		printf("Set environment OMP_NUM_THREADS to 1");
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

	MKL_INT nRhs = 1;

	if (import_with_copy) {
		if(m_dense_values_size > 0)	delete [] m_dense_values;
		if(m_dense_values_fl_size > 0)	delete [] m_dense_values_fl;
	}

	m_dense_values_size = 0;
	m_dense_values_fl_size = 0;

	m_dense_values = NULL;
	m_dense_values_fl = NULL;
}

void DenseSolverMKL::ImportMatrix(SparseMatrix & A) {

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
}

void DenseSolverMKL::ImportMatrix_fl(SparseMatrix & A) {

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
}


void DenseSolverMKL::ImportMatrix_wo_Copy(SparseMatrix & A) {

	USE_FLOAT = false;

	m_rows		= A.rows;
	m_cols		= A.cols;
	m_nnz		= A.nnz;
	m_lda		= A.rows;

	m_dense_values_size = A.dense_values.size();
	m_dense_values_fl_size = 0;

	m_dense_values = &A.dense_values[0];
	
	import_with_copy = false;
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

	if (USE_FLOAT) {
		ssptrf( &U, &m_cols, &m_dense_values_fl[0], &m_ipiv[0] , &info );
	} else {
		dsptrf( &U, &m_cols, &m_dense_values[0], &m_ipiv[0] , &info );
	}

	return info;
}

void DenseSolverMKL::Solve( SEQ_VECTOR <double> & rhs_sol) {

	char U = 'U';
	eslocal info = 0;
	m_nRhs = 1;

	if (USE_FLOAT) {
		tmp_sol_fl.resize(rhs_sol.size());

		for (eslocal i = 0; i < rhs_sol.size(); i++)
			tmp_sol_fl[i] = (float)rhs_sol[i];

		ssptrs( &U, &m_rows, &m_nRhs, &m_dense_values_fl[0], &m_ipiv[0], &tmp_sol_fl[0], &m_rows, &info );	

		for (eslocal i = 0; i < rhs_sol.size(); i++)
			rhs_sol[i] = (double)tmp_sol_fl[i];
	} else {
		dsptrs( &U, &m_rows, &m_nRhs, &m_dense_values[0], &m_ipiv[0], &rhs_sol[0], &m_rows, &info );	
	}
}

void DenseSolverMKL::Solve( SEQ_VECTOR <double> & rhs, SEQ_VECTOR <double> & sol, MKL_INT n_rhs) {

	char U = 'U';
	eslocal info = 0;

	if (USE_FLOAT) {
		sol.resize(rhs.size());
		tmp_sol_fl.resize(rhs.size());

		for (eslocal i = 0; i < rhs.size(); i++)
			tmp_sol_fl[i] = (float)rhs[i];

		ssptrs( &U, &m_rows, &n_rhs, &m_dense_values_fl[0], &m_ipiv[0], &tmp_sol_fl[0], &m_rows, &info );	

		for (eslocal i = 0; i < rhs.size(); i++)
			sol[i] = (double)tmp_sol_fl[i];
	} else {
		sol = rhs;	
		dsptrs( &U, &m_rows, &n_rhs, &m_dense_values[0], &m_ipiv[0], &sol[0], &m_rows, &info );	
	}
}

void DenseSolverMKL::Solve( SEQ_VECTOR <double> & rhs, SEQ_VECTOR <double> & sol, MKL_INT rhs_start_index, MKL_INT sol_start_index) {

	printf("DenseSolverMKL::Solve( SEQ_VECTOR <double> & rhs, SEQ_VECTOR <double> & sol, MKL_INT rhs_start_index, MKL_INT sol_start_index) not implemented yet.\n");
	exit(1);

}