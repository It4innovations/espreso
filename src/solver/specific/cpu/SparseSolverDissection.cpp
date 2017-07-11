
#include "SparseSolverDissection.h"

#include "../../../basis/utilities/utils.h"

#include <Driver/DissectionSolver.hpp>

#ifdef DEBUG
	extern"C" {
	#include <scotch.h>
	}
#endif

using namespace espreso;

SparseSolverDissection::SparseSolverDissection(){

	keep_factors=true;
	initialized = false;
	USE_FLOAT = false;
	import_with_copy = false;

	CSR_I_row_indices_size = 0;
	CSR_J_col_indices_size = 0;
	CSR_V_values_size = 0;
    CSR_V_values_fl_size = 0;
    CSR_I_row_indices = NULL;
    CSR_J_col_indices = NULL;
    CSR_V_values = NULL;
    CSR_V_values_fl = NULL;

	mtype = 2; 			/* Real symmetric positive definite decommatrix */

    // Symbolic refactorization settings
    decomposer = 1; // 0 for SCOTCH, 1 for METIS, 2 for TRIDIAG(Cuthill-McKee)
    eps_pivot = 1.e-2; // pivot threshold
    nb_levels = -1; // number of level of dissection //-1 (automatic)
    scaling = 2; //0;
    kernel_detection_all = true; //false; // is singular?
    min_nodes = 256; // not set
    dslv = NULL;
    diss_verbose = false;
    fp = fopen("/dev/null", "w"); //stdout; // stdout or stderr will cause an error in python benchmarks script
    called = 0;

    MPIrank = 0;
	num_procs = 1;

	rows = 0;
	cols = 0;
	nnz = 0;

	is_whole = false;
	is_sym = true;
	upper_flag = true;

	m_Kplus_size = 0;

#ifdef DEBUG
	// fp = fopen("diss_log_wcp.txt", "w");
	fp = stderr; // stdout or stderr will cause an error in python benchmarks script
	diss_verbose = true;
	int vers, rela, patc;
	SCOTCH_version(&vers, &rela, &patc);

	printf("SCOTCH version: %d.%d.%d\n", vers, rela, patc);
#endif
	error = 0;				/* Initialize error flag */

	/* -------------------------------------------------------------------- */
	/* .. Initialize the internal solver memory pointer. This is only */
	/* necessary for the FIRST call of the PARDISO solver. */
	/* -------------------------------------------------------------------- */

	m_nRhs		 = 1;
	m_factorized = 0;
}

SparseSolverDissection::~SparseSolverDissection() {

		this->Clear();

		if(dslv != NULL) {
			delete dslv;
			dslv = NULL;
		}

		// if(fp != NULL) {
		// 	fclose(fp);
		// 	fp = NULL;
		// }
}

void SparseSolverDissection::Clear() {

	if ( initialized == true )
	{
		MKL_INT nRhs = 1;

		/* -------------------------------------------------------------------- */
		/* .. Termination and release of memory. */
		/* -------------------------------------------------------------------- */

		initialized = false;
	}

	if (import_with_copy) {

		if (CSR_I_row_indices_size > 0) {
			delete [] CSR_I_row_indices;
			CSR_I_row_indices = NULL;
		}
		if (CSR_J_col_indices_size > 0) {
			delete [] CSR_J_col_indices;
			CSR_J_col_indices = NULL;
		}		
		if (CSR_V_values_size > 0) {
			delete [] CSR_V_values;
			CSR_V_values = NULL;
		}	
		//if (CSR_V_values_fl_size > 0)		delete [] CSR_V_values_fl;
	}
	else {
		if (CSR_I_row_indices_size > 0) {
			delete [] CSR_I_row_indices;
			CSR_I_row_indices = NULL;
		}     
		if (CSR_J_col_indices_size > 0)	{
			delete [] CSR_J_col_indices;
			CSR_J_col_indices = NULL;
		}	
	}

	if (USE_FLOAT) {
		if (CSR_V_values_fl_size > 0)		delete [] CSR_V_values_fl;
	}

	CSR_I_row_indices_size = 0;
	CSR_J_col_indices_size = 0;
	CSR_V_values_size      = 0;
	CSR_V_values_fl_size   = 0;

}

void SparseSolverDissection::ImportMatrix(espreso::SparseMatrix & A) {

	USE_FLOAT = false;

	// mkl_set_num_threads();
	is_whole = false;
	is_sym = true;
	upper_flag = true;

	rows	= A.rows;
	cols	= A.cols;
	nnz		= A.nnz;
	m_Kplus_size = A.rows;

	switch (A.mtype) {
	case espreso::MatrixType::REAL_SYMMETRIC_POSITIVE_DEFINITE:
		mtype = 2;
		break;
	case espreso::MatrixType::REAL_SYMMETRIC_INDEFINITE:
		mtype = -2;
		break;
	case espreso::MatrixType::REAL_UNSYMMETRIC:
		mtype = 11;
		is_whole = true;
		is_sym = false;
		break;
	}

	CSR_I_row_indices_size = A.CSR_I_row_indices.size();
	CSR_J_col_indices_size = A.CSR_J_col_indices.size();
	CSR_V_values_size	   = A.CSR_V_values.size();

	CSR_I_row_indices = new MKL_INT[CSR_I_row_indices_size];
	CSR_J_col_indices = new MKL_INT[CSR_J_col_indices_size];
	CSR_V_values	  = new double  [CSR_V_values_size];

	copy(A.CSR_I_row_indices.begin(), A.CSR_I_row_indices.end(), CSR_I_row_indices);
	copy(A.CSR_J_col_indices.begin(), A.CSR_J_col_indices.end(), CSR_J_col_indices);
	copy(A.CSR_V_values     .begin(), A.CSR_V_values     .end(), CSR_V_values);

	// Index base form 1 to 0
	for (int i = 0; i < CSR_I_row_indices_size; ++i)
	{
		CSR_I_row_indices[i] -= 1;
	}

	for (int i = 0; i < CSR_J_col_indices_size; ++i)
	{
		CSR_J_col_indices[i] -= 1;
	}

	import_with_copy = true;

}

void SparseSolverDissection::ImportMatrix_fl(espreso::SparseMatrix & A) {

	printf("Method ImportMatrix_fl is not implemented - float not available in Dissection solver.\n");
	exit(1);

	// USE_FLOAT = true;

	// // sets the num_procs
	// // SetThreaded();
	// // num_procs = 1;
	// int called = 0;

	// dslv_fl = new DissectionSolver<float>(num_procs, diss_verbose, called, fp);

	// // mkl_set_num_threads();
	// is_whole = false;
	// is_sym = true;
	// upper_flag = true;

	// rows	= A.rows;
	// cols	= A.cols;
	// nnz		= A.nnz;
	// m_Kplus_size = A.rows;

	// CSR_I_row_indices_size = A.CSR_I_row_indices.size();
	// CSR_J_col_indices_size = A.CSR_J_col_indices.size();
	// CSR_V_values_size	   = 0;
	// CSR_V_values_fl_size   = A.CSR_V_values.size();

	// CSR_I_row_indices = new MKL_INT[CSR_I_row_indices_size];
	// CSR_J_col_indices = new MKL_INT[CSR_J_col_indices_size];
	// CSR_V_values_fl	  = new float  [CSR_V_values_fl_size];

	// copy(A.CSR_I_row_indices.begin(), A.CSR_I_row_indices.end(), CSR_I_row_indices);
	// copy(A.CSR_J_col_indices.begin(), A.CSR_J_col_indices.end(), CSR_J_col_indices);

	// for (eslocal i = 0; i < CSR_V_values_fl_size; i++)
	// 	CSR_V_values_fl[i] = (float) A.CSR_V_values[i];

	// import_with_copy = true;

}

void SparseSolverDissection::ImportMatrix_wo_Copy_fl(espreso::SparseMatrix & A) {

	printf("Method ImportMatrix_wo_Copy_fl is not implemented - float not available in Dissection solver.\n");
	exit(1);

	// USE_FLOAT = true;

	// rows	= A.rows;
	// cols	= A.cols;
	// nnz		= A.nnz;
	// m_Kplus_size = A.rows;

	// CSR_I_row_indices_size = A.CSR_I_row_indices.size();
	// CSR_J_col_indices_size = A.CSR_J_col_indices.size();
	// CSR_V_values_size	 = 0;
	// CSR_V_values_fl_size = A.CSR_V_values.size();

	// CSR_I_row_indices    = &A.CSR_I_row_indices[0];
	// CSR_J_col_indices    = &A.CSR_J_col_indices[0];

	// CSR_V_values_fl	  = new float  [CSR_V_values_fl_size];
	// for (eslocal i = 0; i < CSR_V_values_fl_size; i++)
	// 	CSR_V_values_fl[i] = (float) A.CSR_V_values[i];

	// import_with_copy = false;
}

void SparseSolverDissection::ImportMatrix_wo_Copy(espreso::SparseMatrix & A) {

	USE_FLOAT = false;

	is_whole = false;
	is_sym = true;
	upper_flag = true;

	rows	= A.rows;
	cols	= A.cols;
	nnz		= A.nnz;
	m_Kplus_size = A.rows;

	switch (A.mtype) {
	case espreso::MatrixType::REAL_SYMMETRIC_POSITIVE_DEFINITE:
		mtype = 2;
		break;
	case espreso::MatrixType::REAL_SYMMETRIC_INDEFINITE:
		mtype = -2;
		break;
	case espreso::MatrixType::REAL_UNSYMMETRIC:
		mtype = 11;
		is_whole = true;
		is_sym = false;
		break;
	}

	CSR_I_row_indices_size = A.CSR_I_row_indices.size();
	CSR_J_col_indices_size = A.CSR_J_col_indices.size();
	CSR_V_values_size	   = A.CSR_V_values.size();

	// CSR_I_row_indices = &A.CSR_I_row_indices[0];
	// CSR_J_col_indices = &A.CSR_J_col_indices[0];
	CSR_V_values	  = &A.CSR_V_values[0];

	CSR_I_row_indices = new MKL_INT[CSR_I_row_indices_size];
	CSR_J_col_indices = new MKL_INT[CSR_J_col_indices_size];

	copy(A.CSR_I_row_indices.begin(), A.CSR_I_row_indices.end(), CSR_I_row_indices);
	copy(A.CSR_J_col_indices.begin(), A.CSR_J_col_indices.end(), CSR_J_col_indices);

	// Index base form 1 to 0
	for (int i = 0; i < CSR_I_row_indices_size; ++i)
	{
		CSR_I_row_indices[i] -= 1;
	}

	for (int i = 0; i < CSR_J_col_indices_size; ++i)
	{
		CSR_J_col_indices[i] -= 1;
	}

	import_with_copy = false;
}

void SparseSolverDissection::ExportMatrix(espreso::SparseMatrix & A) {

	switch (mtype) {
	case 2:
		A.mtype = espreso::MatrixType::REAL_SYMMETRIC_POSITIVE_DEFINITE;
		A.type = 'S';
		break;
	case -2:
		A.mtype = espreso::MatrixType::REAL_SYMMETRIC_INDEFINITE;
		A.type = 'S';
		break;
	case 11:
		A.mtype = espreso::MatrixType::REAL_UNSYMMETRIC;
		A.type = 'G';
		break;
	}

	A.rows = rows;
	A.cols = cols;
	A.nnz = nnz;

	A.CSR_V_values.clear();
	A.CSR_V_values.insert(A.CSR_V_values.end(), &CSR_V_values[0], &CSR_V_values[CSR_V_values_size]);

	A.CSR_I_row_indices.resize(CSR_I_row_indices_size);
	A.CSR_J_col_indices.resize(CSR_J_col_indices_size);

	// Index base from 0 to 1
	for (int i = 0; i < CSR_I_row_indices_size; ++i)
	{
		A.CSR_I_row_indices[i] = CSR_I_row_indices[i] + 1;
	}

	for (int i = 0; i < CSR_J_col_indices_size; ++i)
	{
		A.CSR_J_col_indices[i] = CSR_J_col_indices[i] + 1;
	}
}

void SparseSolverDissection::SetThreaded() {

	/* Numbers of processors, value of OMP_NUM_THREADS */
	num_procs = Esutils::getEnv<int>("SOLVER_NUM_THREADS");
}

int SparseSolverDissection::Factorization(const std::string &str) {

	/* -------------------------------------------------------------------- */
	/* .. Reordering and Symbolic Factorization. This step also allocates */
	/* all memory that is necessary for the factorization. */
	/* -------------------------------------------------------------------- */

	ESINFO(PROGRESS3) << Info::plain() << "f";

	if (USE_FLOAT) {
		printf("Method Factorization for float is not implemented - float not available in Dissection solver.\n");
		exit(1);

		// dslv_fl->SymbolicFact(rows, CSR_I_row_indices, CSR_J_col_indices, is_sym, 
		// 	upper_flag, is_whole, decomposer, nb_levels, min_nodes);
	} else {

		dslv = new DissectionSolver<double>(num_procs, diss_verbose, called, fp);

		dslv->SymbolicFact(rows,(MKL_INT *) CSR_I_row_indices,(MKL_INT *) CSR_J_col_indices, is_sym, upper_flag, is_whole, decomposer, nb_levels, min_nodes);
		
		// dslv->SaveMMMatrix(1, CSR_V_values);
	}

	initialized = true;

	/* -------------------------------------------------------------------- */
	/* .. Numerical factorization. */
	/* -------------------------------------------------------------------- */

	if (USE_FLOAT) {
		printf("Method Factorization for float is not implemented - float not available in Dissection solver.\n");
		exit(1);

	} else {
		dslv->NumericFact(called, (double *) CSR_V_values, scaling, eps_pivot, kernel_detection_all);
	}

	m_factorized = 1;

	//TODO:
	if (USE_FLOAT) {
		tmp_sol_fl1.resize(m_Kplus_size);
		tmp_sol_fl2.resize(m_Kplus_size);
	} else {
		tmp_sol.resize(m_Kplus_size); // - POZOR mozna se musi odkomentovat kvuli alokaci tmp_sol
	}

	return 0;
}

void SparseSolverDissection::Solve( SEQ_VECTOR <double> & rhs_sol) {

	if( USE_FLOAT ) {
		printf("Method Solve for float is not implemented yet - float not available in Dissection solver.\n");
		exit(1);

		// for (eslocal i = 0; i < m_Kplus_size; i++)
		// 	tmp_sol_fl1[i] = (float)rhs_sol[i];
	}

	if (!initialized) {
		std::stringstream ss;
		ss << "Solve -> rank: " << environment->MPIrank;
		Factorization(ss.str());
	}

	MKL_INT n_rhs = 1;

	/* -------------------------------------------------------------------- */
	/* .. Back substitution and iterative refinement. */
	/* -------------------------------------------------------------------- */
	
	if (USE_FLOAT) {
		printf("Method Solve for float is not implemented - float not available in Dissection solver.\n");
		exit(1);
		
	} else {

		bool projection = true;
		bool is_trans = false;
		bool is_scaling = true;
		dslv->SolveSingle(&rhs_sol.front(), projection, is_trans, is_scaling);

		// projection -> až se bude přebírat jádro
		// scaling 0 - no scaling
		// 		   1 - sedlobodove systemy
		// 		   2 - SP semi D
	}

	if (!keep_factors) {
		/* -------------------------------------------------------------------- */
		/* .. Termination and release of memory. */
		/* -------------------------------------------------------------------- */
		dslv->NumericFree();
		
		initialized = false;
	}

	if( USE_FLOAT ) {
		printf("Method Solve for float is not implemented - float not available in Dissection solver.\n");
		exit(1);
		// for (eslocal i = 0; i < m_Kplus_size; i++)
		// 	rhs_sol[i] = (double)tmp_sol_fl1[i];
	}
}

void SparseSolverDissection::Solve( SEQ_VECTOR <double> & rhs, SEQ_VECTOR <double> & sol, MKL_INT n_rhs) {

	SEQ_VECTOR <float> tmp_in, tmp_out;

	if( USE_FLOAT ) {
		printf("Method Solve for float is not implemented yet - float not available in Dissection solver.\n");
		exit(1);

		// tmp_in.resize(n_rhs * m_Kplus_size);
		// tmp_out.resize(n_rhs * m_Kplus_size);

		// for (eslocal i = 0; i < n_rhs * m_Kplus_size; i++)
		// 	tmp_in[i] = (float)rhs[i];
	}

	if (!initialized) {
		std::stringstream ss;
		ss << "Solve -> rank: " << environment->MPIrank;
		Factorization(ss.str());
	}

	/* -------------------------------------------------------------------- */
	/* .. Back substitution and iterative refinement. */
	/* -------------------------------------------------------------------- */

	if (USE_FLOAT) {
		printf("Method Solve for float is not implemented yet - float not available in Dissection solver.\n");
		exit(1);

	} else {
		bool projection = true;
		bool is_trans = false;
		bool is_scaling = true;

		sol = rhs;

		// Alternative multiple RHS solve
//		for(eslocal i=0; i< n_rhs; i++) {
//			dslv->SolveSingle(&sol[i*rows], projection, is_trans, is_scaling);
//		}

		dslv->SolveMulti(&sol.front(), n_rhs, projection, is_trans, is_scaling);
	}

	if (!keep_factors) {
		/* -------------------------------------------------------------------- */
		/* .. Termination and release of memory. */
		/* -------------------------------------------------------------------- */
		dslv->NumericFree();
		
		initialized = false;
	}

	if( USE_FLOAT ) {
		printf("Method Solve for float is not implemented yet - float not available in Dissection solver.\n");
		exit(1);

		// for (eslocal i = 0; i < n_rhs * m_Kplus_size; i++)
		// 	sol[i] = (double)tmp_out[i];
	}
}

void SparseSolverDissection::Solve( SEQ_VECTOR <double> & rhs, SEQ_VECTOR <double> & sol, MKL_INT rhs_start_index, MKL_INT sol_start_index) {

	if( USE_FLOAT ) {
		printf("Method Solve for float is not implemented yet.\n");
		exit(1);

		// for (eslocal i = 0; i < m_Kplus_size; i++)
		// 	tmp_sol_fl1[i] = (float)rhs[rhs_start_index+ i];
	}

	if (!initialized) {
		std::stringstream ss;
		ss << "Solve -> rank: " << environment->MPIrank;
		Factorization(ss.str());
	}

	/* -------------------------------------------------------------------- */
	/* .. Back substitution and iterative refinement. */
	/* -------------------------------------------------------------------- */

	if (USE_FLOAT) {
		printf("Method Solve for float is not implemented yet - float not available in Dissection solver.\n");
		exit(1);

	} else {

		bool projection = true;
		bool is_trans = false;
		bool is_scaling = true;

		copy(&rhs[rhs_start_index], &rhs[rhs_start_index + rows], &sol[sol_start_index]);

		dslv->SolveSingle(&sol[sol_start_index], projection, is_trans, is_scaling);
	}

	if (!keep_factors) {
		/* -------------------------------------------------------------------- */
		/* .. Termination and release of memory. */
		/* -------------------------------------------------------------------- */
		dslv->NumericFree();
		
		initialized = false;
	}

	if( USE_FLOAT ) {
		printf("Method Solve for float is not implemented yet - float not available in Dissection solver.\n");
		exit(1);
		// for (eslocal i = 0; i < m_Kplus_size; i++)
		// 	sol[i + sol_start_index] = (double)tmp_sol_fl2[i];
	}
}

void SparseSolverDissection::SolveMat_Sparse( espreso::SparseMatrix & A) {
	SolveMat_Sparse(A, A);
};

void SparseSolverDissection::SolveMat_Sparse( espreso::SparseMatrix & A_in, espreso::SparseMatrix & B_out) {
	SolveMat_Sparse(A_in, B_out, 'T');
};

void SparseSolverDissection::SolveMat_Sparse( espreso::SparseMatrix & A_in, espreso::SparseMatrix & B_out, char T_for_input_matrix_is_transposed_N_input_matrix_is_NOT_transposed ) {

	if (!initialized) {
		std::stringstream ss;
		ss << "Solve -> rank: " << environment->MPIrank;
		Factorization(ss.str());
	}

	bool keep_factors_tmp = keep_factors;
	keep_factors          = true;

	char trans = T_for_input_matrix_is_transposed_N_input_matrix_is_NOT_transposed;

	espreso::SparseMatrix tmpM;
	if (trans == 'T')
		A_in.MatTranspose(tmpM);
	else
		tmpM = A_in;


	SEQ_VECTOR<double> rhs;
	SEQ_VECTOR<double> sol;

	for (size_t row = 1, offset = 0; row < tmpM.CSR_I_row_indices.size(); row++) {
		if (tmpM.CSR_I_row_indices[row - 1] < tmpM.CSR_I_row_indices[row]) {
			offset = rhs.size();
			rhs.resize(rhs.size() + tmpM.cols);
			for (eslocal col = tmpM.CSR_I_row_indices[row - 1]; col < tmpM.CSR_I_row_indices[row]; col++) {
				rhs[offset + tmpM.CSR_J_col_indices[col - 1] - 1] = tmpM.CSR_V_values[col - 1];
			}
		}
	}

	sol.resize(rhs.size());
	if (tmpM.cols !=0) {
		Solve(rhs, sol, rhs.size() / tmpM.cols);
	}
	
	for (size_t row = 1, offset = 0; row < tmpM.CSR_I_row_indices.size(); row++) {
		if (tmpM.CSR_I_row_indices[row - 1] < tmpM.CSR_I_row_indices[row]) {
			for (eslocal col = 0; col < tmpM.cols; col++, offset++){
				if (sol[offset] != 0.0) {
					tmpM.I_row_indices.push_back(row);
					tmpM.J_col_indices.push_back(col + 1);
					tmpM.V_values.push_back(sol[offset]);
				}
			}
		}
	}

	tmpM.nnz = tmpM.V_values.size();
	tmpM.ConvertToCSR(1);
	tmpM.MatTranspose(B_out);

	tmpM.Clear();

	keep_factors = keep_factors_tmp;
	if (!keep_factors) {
		/* -------------------------------------------------------------------- */
		/* .. Termination and release of memory. */
		/* -------------------------------------------------------------------- */
//		phase = -1;			/* Release internal memory. */
//		MKL_INT nRhs = 1;
//		double ddum;			/* Double dummy */
//		MKL_INT idum;			/* Integer dummy. */
//		PARDISO (pt, &maxfct, &mnum, &mtype, &phase,
//				&rows, &ddum, CSR_I_row_indices, CSR_J_col_indices, &idum, &nRhs,
//				iparm, &msglvl, &ddum, &ddum, &error);
		initialized = false;
	}
}


void SparseSolverDissection::SolveMat_Dense( espreso::SparseMatrix & A ) {
	SolveMat_Dense(A, A);
}

void SparseSolverDissection::SolveMat_Dense( espreso::SparseMatrix & A_in, espreso::SparseMatrix & B_out ) {

	SEQ_VECTOR<double> rhs;
	SEQ_VECTOR<double> sol;

	MKL_INT job[8];
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

	MKL_INT lda  = m;
	MKL_INT info;


	// Convert input matrix (RHS) to dense format

	//void mkl_ddnscsr (
	//	int *job,
	//	int *m, int *n,
	//	double *Adns, int *lda,
	//	double *Acsr, int *AJ, int *AI,
	//	int *info);

	mkl_ddnscsr (
		job,
		&m, &n,
		&rhs[0], &lda,
		&A_in.CSR_V_values[0], &A_in.CSR_J_col_indices[0], &A_in.CSR_I_row_indices[0],
		&info);

	// Solve with multiple right hand sides
	Solve(rhs,sol,nRhs);
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
	MKL_INT nnzmax = B_out.CSR_I_row_indices[m];//-1;

	B_out.CSR_J_col_indices.resize(nnzmax);
	B_out.CSR_V_values.resize(nnzmax);

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

	sol.clear();

	// Setup parameters for output matrix
	B_out.cols	= A_in.cols;
	B_out.rows	= A_in.rows;
	B_out.nnz	= B_out.CSR_V_values.size();
	B_out.type	= 'G';

}

//Obsolete - to be removed
void SparseSolverDissection::SolveMatF( espreso::SparseMatrix & A_in, espreso::SparseMatrix & B_out, bool isThreaded ) {

	printf("Method SolveMatF is not implemented yet.\n");
	exit(1);	
}

void SparseSolverDissection::Create_SC( espreso::SparseMatrix & SC_out, MKL_INT sc_size, bool isThreaded ) {

	SparseMatrix K_modif;
	SparseMatrix K_rr;
	SparseMatrix K_rs;
	SparseMatrix K_sr;
	SparseMatrix KsrInvKrrKrs;

	ExportMatrix(K_modif); // memory inefficient

	eslocal i_start = 0;
	eslocal nonsing_size = K_modif.rows - sc_size - i_start;
	eslocal j_start = nonsing_size;

	K_rs.getSubBlockmatrix_rs(K_modif,K_rs,i_start, nonsing_size,j_start,sc_size);

	if (mtype == 2){ // symmetric system
		K_rr.mtype = espreso::MatrixType::REAL_SYMMETRIC_POSITIVE_DEFINITE;
        K_rs.MatTranspose(K_sr);
	} else { // unsymmetric system
		K_sr.getSubBlockmatrix_rs(K_modif,K_sr,j_start,sc_size,i_start, nonsing_size);
	}

	SC_out.getSubDiagBlockmatrix(K_modif,SC_out,nonsing_size,sc_size);

	// K_rs =    inv(K_rr) * K_rs          (classical Dirichlet precond. assembled by own - not via PardisoSC routine)
	K_rr.getSubDiagBlockmatrix(K_modif,K_rr,i_start, nonsing_size);

	SparseSolverDissection K_rr_solver;
	K_rr_solver.ImportMatrix_wo_Copy(K_rr);
	K_rr_solver.SolveMat_Dense(K_rs);

	KsrInvKrrKrs.MatMat(K_sr,'N',K_rs);
	SC_out.MatAddInPlace(KsrInvKrrKrs,'N',-1);

    // TODO - changed for AD-system
    if ( mtype == 11){
      SC_out.type = 'G';
      SC_out.MatTranspose();
    }
    else
    {
      SC_out.RemoveLower();
      SC_out.type = 'S';
    }

	if (msglvl) {
		ESINFO(EXHAUSTIVE) << SC_out.SpyText();
	}

//    if (generate_symmetric_sc_1_generate_general_sc_0 == 1) {
//    	SC_out.RemoveLower();
//    }
}

void SparseSolverDissection::Create_SC_w_Mat( espreso::SparseMatrix & K_in, espreso::SparseMatrix & B_in, espreso::SparseMatrix & SC_out,
								    bool isThreaded, MKL_INT generate_symmetric_sc_1_generate_general_sc_0 ) {

	// |  K_in      B_in |
	// | (B_in)t     0   |

	SparseMatrix K_sc1;
	SparseMatrix K_b_tmp;

	SparseSolverDissection K_sc1_solver;
	K_sc1_solver.ImportMatrix_wo_Copy(K_in);
	if(isThreaded)
		K_sc1_solver.SetThreaded();
	K_sc1_solver.SolveMat_Dense(B_in, K_b_tmp);

	SC_out.MatMat(B_in,'T',K_b_tmp);

	//SC_out.MatScale(-1.0);

	SC_out.type = 'G';

	if (msglvl) {
		ESINFO(EXHAUSTIVE) << SC_out.SpyText();
	}

	SC_out.ConvertCSRToDense(1);

	if (generate_symmetric_sc_1_generate_general_sc_0 == 1) {
		//SC_out.RemoveLower();
		SC_out.RemoveLowerDense();
		SC_out.type = 'S';
	}
}

void SparseSolverDissection::Create_non_sym_SC_w_Mat( espreso::SparseMatrix & K_in, espreso::SparseMatrix & B1_in, espreso::SparseMatrix & B0_in, espreso::SparseMatrix & SC_out, bool isThreaded, MKL_INT generate_symmetric_sc_1_generate_general_sc_0 ) {

	// |  K_in      B1_in |
	// | (B0_in)t     0   |

	ESINFO(PROGRESS1) << K_in.SpyText();
	ESINFO(PROGRESS1) << B1_in.SpyText();
	ESINFO(PROGRESS1) << B0_in.SpyText();

	SparseMatrix K_sc1;
	SparseMatrix K_b_tmp;

	SparseSolverDissection K_sc1_solver;
	K_sc1_solver.ImportMatrix_wo_Copy(K_in);
	if(isThreaded)
		K_sc1_solver.SetThreaded();
	K_sc1_solver.SolveMat_Dense(B1_in, K_b_tmp);

	SC_out.MatMat(B0_in,'T',K_b_tmp);

	//SC_out.MatScale(-1.0);

	SC_out.type = 'G';

	if (msglvl) {
		ESINFO(EXHAUSTIVE) << SC_out.SpyText();
	}

	SC_out.ConvertCSRToDense(0);

	if (generate_symmetric_sc_1_generate_general_sc_0 == 1) {
		SC_out.RemoveLowerDense();
		SC_out.RemoveLower();
		SC_out.type = 'S';
	}
}

void SparseSolverDissection::GetKernel(SparseMatrix &R) {

	// After factorization, with singular matrix
	if (!initialized) {
		std::stringstream ss;
		ss << "Get Kernel -> rank: " << environment->MPIrank;
		Factorization(ss.str());
	}

	eslocal kern_dim = dslv->kern_dimension();
	R.dense_values.resize(kern_dim * rows);

	dslv->GetKernelVectors(&R.dense_values.front());

	R.type = 'G';
	R.nnz  = kern_dim * cols;

	R.cols = kern_dim;
	R.rows = rows;




#if 1

	SEQ_VECTOR <double> kern_vec;
	kern_vec = R.dense_values;

	// Test the obtained kernel vectors
	// Test 1
	// A*r = 0
	ESINFO(PROGRESS1) << "Dissection kernel test 1: 'A*r = 0'";
	SEQ_VECTOR <double> r(kern_dim * rows, -1.0);
	SparseMatrix K_sing;
	ExportMatrix(K_sing);

	for(eslocal i = 0; i < kern_dim; i++) {
		K_sing.MatVec(kern_vec, r, 'N', i*rows, i*rows);
		double norm_r = 0.0;
		for(int j = i*rows; j < i*rows+rows; j++) {
			norm_r += r[j] * r[j];
		}
		norm_r = sqrt(norm_r);
		ESINFO(PROGRESS1) << "||r" << i << "|| = " << norm_r;
	}

	// Test 2
	// ||K * R|| / (||K|| * ||R||)
	SparseMatrix Rr;
	Rr.dense_values = kern_vec;
	Rr.cols = kern_dim;
	Rr.rows = kern_vec.size() / kern_dim;

	double norm_K_R = K_sing.getNorm_K_R(K_sing,Rr,'N');

	double norm_R = 0.0;
	for(int i = 0; i < Rr.dense_values.size(); i++) {
		norm_R += Rr.dense_values[i] * Rr.dense_values[i];
	}
	norm_R = sqrt(norm_R);

	K_sing.ConvertCSRToDense(0);
	double norm_K = 0.0;
	for(int i = 0; i < K_sing.dense_values.size(); i++) {
		norm_K += K_sing.dense_values[i] * K_sing.dense_values[i];
	}
	norm_K = sqrt(norm_K);

	double test2 = norm_K_R / (norm_K * norm_R);
	ESINFO(PROGRESS1) << "\nDissection kernel test 2: '||K * R|| / (||K|| * ||R||)'";
	ESINFO(PROGRESS1) << test2;

	// Test 3
	// R^T * R = identity
	SparseMatrix I;
	Rr.ConvertDenseToCSR(0);
	I.MatMat(Rr, 'T', Rr);
	ESINFO(PROGRESS1) << "\nDissection kernel test 3: 'R^T * R = I'";
	ESINFO(PROGRESS1) << I;

	// Test 4
	// ||K * R|| / max(diag(K))
	double max_diag_K = K_sing.getDiagonalMaximum();
	double test4 = norm_K_R / max_diag_K;
	ESINFO(PROGRESS1) << "Dissection kernel test 4: '||K * R|| / max(diag(K))'";
	ESINFO(PROGRESS1) << test4 << "\n";
#endif


}

void SparseSolverDissection::GetKernelVectors(SEQ_VECTOR <double> & kern_vec, eslocal & kern_dim) {

	// After factorization, with singular matrix
	if (!initialized) {
		std::stringstream ss;
		ss << "Get Kernel -> rank: " << environment->MPIrank;
		Factorization(ss.str());
	}

	kern_dim = dslv->kern_dimension();
	kern_vec.resize(kern_dim * rows);
	dslv->GetKernelVectors(&kern_vec.front());

#if 1
	// Test the obtained kernel vectors
	// Test 1
	// A*r = 0
	ESINFO(PROGRESS1) << "Dissection kernel test 1: 'A*r = 0'";
	SEQ_VECTOR <double> r(kern_dim * rows, -1.0);
	SparseMatrix K_sing;
	ExportMatrix(K_sing);

	for(eslocal i = 0; i < kern_dim; i++) {
		K_sing.MatVec(kern_vec, r, 'N', i*rows, i*rows);
		double norm_r = 0.0;
		for(int j = i*rows; j < i*rows+rows; j++) {
			norm_r += r[j] * r[j];
		}
		norm_r = sqrt(norm_r);
		ESINFO(PROGRESS1) << "||r" << i << "|| = " << norm_r;
	}

	// Test 2
	// ||K * R|| / (||K|| * ||R||)
	SparseMatrix R;
	R.dense_values = kern_vec;
	R.cols = kern_dim;
	R.rows = kern_vec.size() / kern_dim;

	double norm_K_R = K_sing.getNorm_K_R(K_sing,R,'N');

	double norm_R = 0.0;
	for(int i = 0; i < R.dense_values.size(); i++) {
		norm_R += R.dense_values[i] * R.dense_values[i];
	}
	norm_R = sqrt(norm_R);

	K_sing.ConvertCSRToDense(0);
	double norm_K = 0.0;
	for(int i = 0; i < K_sing.dense_values.size(); i++) {
		norm_K += K_sing.dense_values[i] * K_sing.dense_values[i];
	}
	norm_K = sqrt(norm_K);

	double test2 = norm_K_R / (norm_K * norm_R);
	ESINFO(PROGRESS1) << "\nDissection kernel test 2: '||K * R|| / (||K|| * ||R||)'";
	ESINFO(PROGRESS1) << test2;

	// Test 3
	// R^T * R = identity
	SparseMatrix I;
	R.ConvertDenseToCSR(0);
	I.MatMat(R, 'T', R);
	ESINFO(PROGRESS1) << "\nDissection kernel test 3: 'R^T * R = I'";
	ESINFO(PROGRESS1) << I;

	// Test 4
	// ||K * R|| / max(diag(K))
	double max_diag_K = K_sing.getDiagonalMaximum();
	double test4 = norm_K_R / max_diag_K;
	ESINFO(PROGRESS1) << "Dissection kernel test 4: '||K * R|| / max(diag(K))'";
	ESINFO(PROGRESS1) << test4 << "\n";
#endif
}


#include <stdio.h>
#include "mkl_rci.h"
#include "mkl_blas.h"
#include "mkl_spblas.h"
#include "mkl_service.h"

void SparseSolverDissection::SolveCG(espreso::SparseMatrix & A_in, SEQ_VECTOR <double> & rhs_sol) {

	MKL_INT size = A_in.rows;
	SEQ_VECTOR <double> sol (size, 0);

	SolveCG(A_in, rhs_sol, sol);

	rhs_sol = sol;
}

void SparseSolverDissection::SolveCG(espreso::SparseMatrix & A_in, SEQ_VECTOR <double> & rhs_in, SEQ_VECTOR <double> & sol) {
	SEQ_VECTOR<double> init;
	SolveCG(A_in, rhs_in, sol, init);
}

void SparseSolverDissection::SolveCG(espreso::SparseMatrix & A_in, SEQ_VECTOR <double> & rhs_in, SEQ_VECTOR <double> & sol, SEQ_VECTOR <double> & initial_guess) {

	printf("Method SolveCG is not implemented yet.\n");
	exit(1);
}
