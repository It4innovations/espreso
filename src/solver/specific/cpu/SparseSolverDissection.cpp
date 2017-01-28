
#include "SparseSolverDissection.h"

#include "../../../basis/utilities/utils.h"

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
    nb_levels = -1; //8; // number of level of dissection //-1 (automatic)
    scaling = 0;
    kernel_detection_all = false; // is singular?
    min_nodes = 256; // not set
    dslv = NULL;
    // dslv_fl = NULL;
    diss_verbose = false;
    fp = NULL;
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

void SparseSolverDissection::ImportMatrix(SparseMatrix & A) {

	USE_FLOAT = false;

	// num_procs = 1;
	called = 0;

	#ifdef DEBUG
		// fp = fopen("diss_log_cp.txt", "w");
		fp = stderr;
	#endif

	dslv = new DissectionSolver<double>(num_procs, diss_verbose, called, fp);

	// mkl_set_num_threads();
	is_whole = false;
	is_sym = true;
	upper_flag = true;

	rows	= A.rows;
	cols	= A.cols;
	nnz		= A.nnz;
	m_Kplus_size = A.rows;

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

void SparseSolverDissection::ImportMatrix_fl(SparseMatrix & A) {

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

void SparseSolverDissection::ImportMatrix_wo_Copy_fl(SparseMatrix & A) {

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

void SparseSolverDissection::ImportMatrix_wo_Copy(SparseMatrix & A) {

	USE_FLOAT = false;

	// num_procs = 1;
	called = 0;

	#ifdef DEBUG
		// fp = fopen("diss_log_wcp.txt", "w");
		fp = stderr;
	#endif

	dslv = new DissectionSolver<double>(num_procs, diss_verbose, called, fp);

	is_whole = false;
	is_sym = true;
	upper_flag = true;

	rows	= A.rows;
	cols	= A.cols;
	nnz		= A.nnz;
	m_Kplus_size = A.rows;

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

void SparseSolverDissection::SetThreaded() {

	/* Numbers of processors, value of OMP_NUM_THREADS */
	num_procs = Esutils::getEnv<int>("SOLVER_NUM_THREADS");
}

int SparseSolverDissection::Factorization(const std::string &str) {

	/* -------------------------------------------------------------------- */
	/* .. Reordering and Symbolic Factorization. This step also allocates */
	/* all memory that is necessary for the factorization. */
	/* -------------------------------------------------------------------- */

	ESINFO(PROGRESS2) << Info::plain() << "f";

	if (USE_FLOAT) {
		printf("Method Factorization for float is not implemented - float not available in Dissection solver.\n");
		exit(1);

		// dslv_fl->SymbolicFact(rows, CSR_I_row_indices, CSR_J_col_indices, is_sym, 
		// 	upper_flag, is_whole, decomposer, nb_levels, min_nodes);
	} else {

		dslv->SymbolicFact(rows,(MKL_INT *) CSR_I_row_indices,(MKL_INT *) CSR_J_col_indices, is_sym, upper_flag, is_whole, decomposer, nb_levels, min_nodes);
		
		// dslv->SaveMMMatrix(1, CSR_V_values);
	}

	if (error != 0)
	{
    return error;
		SparseMatrix s;
		s.rows = rows;
		s.cols = cols;
		s.type = 'S';
		s.nnz = nnz;
		s.CSR_I_row_indices = std::vector<eslocal>(CSR_I_row_indices, CSR_I_row_indices + CSR_I_row_indices_size);
		s.CSR_J_col_indices = std::vector<eslocal>(CSR_J_col_indices, CSR_J_col_indices + CSR_J_col_indices_size);
		s.CSR_V_values = std::vector<double>(CSR_V_values, CSR_V_values + CSR_V_values_size);

		std::ofstream osK(Logging::prepareFile("ERROR").c_str());
		osK << s;
		osK.close();

		ESINFO(ERROR) << error << " during symbolic factorization";
		exit (EXIT_FAILURE);
	} else {
		initialized = true;
	}

	/* -------------------------------------------------------------------- */
	/* .. Numerical factorization. */
	/* -------------------------------------------------------------------- */

	if (USE_FLOAT) {
		printf("Method Factorization for float is not implemented - float not available in Dissection solver.\n");
		exit(1);

	} else {
		dslv->NumericFact(called, (double *) CSR_V_values, scaling, eps_pivot, kernel_detection_all);
	}

	if (error != 0)
	{
		return error;
		SparseMatrix s;
		s.rows = rows;
		s.cols = cols;
		s.type = 'S';
		s.nnz = nnz;
		s.CSR_I_row_indices = std::vector<eslocal>(CSR_I_row_indices, CSR_I_row_indices + CSR_I_row_indices_size);
		s.CSR_J_col_indices = std::vector<eslocal>(CSR_J_col_indices, CSR_J_col_indices + CSR_J_col_indices_size);
		s.CSR_V_values = std::vector<double>(CSR_V_values, CSR_V_values + CSR_V_values_size);

		std::ofstream osK(Logging::prepareFile("ERROR").c_str());
		osK << s;
		osK.close();

		ESINFO(ERROR) << error << " during numerical factorization";
		exit (EXIT_FAILURE);
	} else {
		m_factorized = 1;
	}
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

		bool projection = false;
		bool is_trans = false;
		bool is_scaling = true;
		dslv->SolveSingle(&rhs_sol.front(), projection, is_trans, is_scaling);

		// projection -> až se bude přebírat jádro
		// scaling 0 - no scaling
		// 		   1 - sedlobodove systemy
		// 		   2 - SP semi D
	}

	if (error != 0)
	{
		SparseMatrix s;
		s.rows = rows;
		s.cols = cols;
		s.type = 'S';
		s.nnz = nnz;
		s.CSR_I_row_indices = std::vector<eslocal>(CSR_I_row_indices, CSR_I_row_indices + CSR_I_row_indices_size);
		s.CSR_J_col_indices = std::vector<eslocal>(CSR_J_col_indices, CSR_J_col_indices + CSR_J_col_indices_size);
		s.CSR_V_values = std::vector<double>(CSR_V_values, CSR_V_values + CSR_V_values_size);

		std::ofstream osK(Logging::prepareFile("ERROR").c_str());
		osK << s;
		osK.close();

		ESINFO(ERROR) << "ERROR during solution: " << error;
		exit (3);
	}

	if (!keep_factors) {
		/* -------------------------------------------------------------------- */
		/* .. Termination and release of memory. */
		/* -------------------------------------------------------------------- */
		MKL_INT nRhs = 1;
		
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
		bool projection = false;
		bool is_trans = false;
		bool is_scaling = true;

		sol = rhs;

		dslv->SolveMulti(&sol.front(), n_rhs, projection, is_trans, is_scaling);
	}

	if (error != 0)
	{
		SparseMatrix s;
		s.rows = rows;
		s.cols = cols;
		s.type = 'S';
		s.nnz = nnz;
		s.CSR_I_row_indices = std::vector<eslocal>(CSR_I_row_indices, CSR_I_row_indices + CSR_I_row_indices_size);
		s.CSR_J_col_indices = std::vector<eslocal>(CSR_J_col_indices, CSR_J_col_indices + CSR_J_col_indices_size);
		s.CSR_V_values = std::vector<double>(CSR_V_values, CSR_V_values + CSR_V_values_size);

		std::ofstream osK(Logging::prepareFile("ERROR").c_str());
		osK << s;
		osK.close();

		ESINFO(ERROR) << "ERROR during solution: " << error;
		exit (3);
	}

	if (!keep_factors) {
		/* -------------------------------------------------------------------- */
		/* .. Termination and release of memory. */
		/* -------------------------------------------------------------------- */
		MKL_INT nRhs = 1;
		
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

		bool projection = false;
		bool is_trans = false;
		bool is_scaling = true;
		dslv->SolveSingle(&rhs[rhs_start_index], projection, is_trans, is_scaling);

	}

	if (error != 0)
	{
		SparseMatrix s;
		s.rows = rows;
		s.cols = cols;
		s.type = 'S';
		s.nnz = nnz;
		s.CSR_I_row_indices = std::vector<eslocal>(CSR_I_row_indices, CSR_I_row_indices + CSR_I_row_indices_size);
		s.CSR_J_col_indices = std::vector<eslocal>(CSR_J_col_indices, CSR_J_col_indices + CSR_J_col_indices_size);
		s.CSR_V_values = std::vector<double>(CSR_V_values, CSR_V_values + CSR_V_values_size);

		std::ofstream osK(Logging::prepareFile("ERROR").c_str());
		osK << s;
		osK.close();

		ESINFO(ERROR) << "ERROR during solution: " << error;
		exit (3);
	}


	if (!keep_factors) {
		/* -------------------------------------------------------------------- */
		/* .. Termination and release of memory. */
		/* -------------------------------------------------------------------- */
		MKL_INT nRhs = 1;
		
		initialized = false;
	}


	if( USE_FLOAT ) {
		printf("Method Solve for float is not implemented yet - float not available in Dissection solver.\n");
		exit(1);
		// for (eslocal i = 0; i < m_Kplus_size; i++)
		// 	sol[i + sol_start_index] = (double)tmp_sol_fl2[i];
	}

}

void SparseSolverDissection::SolveMat_Sparse( SparseMatrix & A) {
	SolveMat_Sparse(A, A);
};

void SparseSolverDissection::SolveMat_Sparse( SparseMatrix & A_in, SparseMatrix & B_out) {
	SolveMat_Sparse(A_in, B_out, 'T');
};

void SparseSolverDissection::SolveMat_Sparse( SparseMatrix & A_in, SparseMatrix & B_out, char T_for_input_matrix_is_transposed_N_input_matrix_is_NOT_transposed ) {

	if (!initialized) {
		std::stringstream ss;
		ss << "Solve -> rank: " << environment->MPIrank;
		Factorization(ss.str());
	}

	bool keep_factors_tmp = keep_factors;
	keep_factors          = true;

	char trans = T_for_input_matrix_is_transposed_N_input_matrix_is_NOT_transposed;

	SparseMatrix tmpM;
	if (trans == 'T')
		A_in.MatTranspose(tmpM);
	else
		tmpM = A_in;


	SEQ_VECTOR<double> rhs;
	SEQ_VECTOR<double> sol;

	rhs.resize(tmpM.cols);
	sol.resize(tmpM.cols);

	// main loop over rows
	MKL_INT col = 0;
	MKL_INT n_nnz = 0;
	for (MKL_INT row = 1; row < tmpM.CSR_I_row_indices.size(); row++) {
		MKL_INT row_size = tmpM.CSR_I_row_indices[row] - tmpM.CSR_I_row_indices[row-1];
		if (row_size > 0) {
			for (MKL_INT c = 0; c < row_size; c++) { // loop over selected row
				rhs[ tmpM.CSR_J_col_indices[col] - 1] = tmpM.CSR_V_values [col];
				col++;
			}
			MKL_INT nRhs_l = 1;
			//m_error = dss_solve_real (m_handle, m_opt, &rhs[0], nRhs_l, &sol[0]);
			Solve(rhs, sol, nRhs_l);

			for (MKL_INT s = 0; s < sol.size(); s++){
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

	keep_factors = keep_factors_tmp;
	if (!keep_factors) {
		/* -------------------------------------------------------------------- */
		/* .. Termination and release of memory. */
		/* -------------------------------------------------------------------- */
		MKL_INT nRhs = 1;
		
		initialized = false;
	}
}


void SparseSolverDissection::SolveMat_Dense( SparseMatrix & A ) {
	SolveMat_Dense(A, A);
}

void SparseSolverDissection::SolveMat_Dense( SparseMatrix & A_in, SparseMatrix & B_out ) {

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
void SparseSolverDissection::SolveMatF( SparseMatrix & A_in, SparseMatrix & B_out, bool isThreaded ) {

	printf("Method SolveMatF is not implemented yet.\n");
	exit(1);	
}

void SparseSolverDissection::Create_SC( SparseMatrix & SC_out, MKL_INT sc_size, bool isThreaded ) {

	printf("Method Create_SC is not implemented yet.\n");
	exit(1);
}


void SparseSolverDissection::Create_SC_w_Mat( SparseMatrix & K_in, SparseMatrix & B_in, SparseMatrix & SC_out,
								    bool isThreaded, MKL_INT generate_symmetric_sc_1_generate_general_sc_0 ) {

	printf("Method Create_SC_w_Mat is not implemented yet.\n");
	exit(1);
}

void SparseSolverDissection::Create_non_sym_SC_w_Mat( SparseMatrix & K_in, SparseMatrix & B1_in, SparseMatrix & B0_in, SparseMatrix & SC_out, bool isThreaded, MKL_INT generate_symmetric_sc_1_generate_general_sc_0 ) {

	printf("Method Create_non_sym_SC_w_Mat is not implemented yet.\n");
	exit(1);
}


#include <stdio.h>
#include "mkl_rci.h"
#include "mkl_blas.h"
#include "mkl_spblas.h"
#include "mkl_service.h"

void SparseSolverDissection::SolveCG(SparseMatrix & A_in, SEQ_VECTOR <double> & rhs_sol) {

	MKL_INT size = A_in.rows;
	SEQ_VECTOR <double> sol (size, 0);

	SolveCG(A_in, rhs_sol, sol);

	rhs_sol = sol;
}

void SparseSolverDissection::SolveCG(SparseMatrix & A_in, SEQ_VECTOR <double> & rhs_in, SEQ_VECTOR <double> & sol) {
	SEQ_VECTOR<double> init;
	SolveCG(A_in, rhs_in, sol, init);
}

void SparseSolverDissection::SolveCG(SparseMatrix & A_in, SEQ_VECTOR <double> & rhs_in, SEQ_VECTOR <double> & sol, SEQ_VECTOR <double> & initial_guess) {

	printf("Method SolveCG is not implemented yet.\n");
	exit(1);
}
