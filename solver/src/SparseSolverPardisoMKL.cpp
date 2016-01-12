#include "SparseSolver.h"
#include "mkl_pardiso.h"












SparseSolver::SparseSolver(){

	keep_factors=true;
	initialized = false;

	CSR_I_row_indices_size = 0;
	CSR_J_col_indices_size = 0;
	CSR_V_values_size = 0;


	mtype = 2; 			/* Real symmetric positive definite matrix */

	/* -------------------------------------------------------------------- */
	/* .. Setup Pardiso control parameters. */
	/* -------------------------------------------------------------------- */

	for (int i = 0; i < 64; i++)
	{
		iparm[i] = 0;
	}
	iparm[0] = 1;			/* No solver default */
	iparm[1] = 2;			/* Fill-in reordering from METIS */

//	/* Numbers of processors, value of OMP_NUM_THREADS */
//	int num_procs;
//	char * var = getenv("OMP_NUM_THREADS");
//    if(var != NULL)
//    	sscanf( var, "%d", &num_procs );
//	else {
//    	printf("Set environment OMP_NUM_THREADS to 1");
//        exit(1);
//	}
//    iparm[2]  = num_procs;

	iparm[2] = 1; 			//by default the solver runs with single thread

	iparm[3] = 0;//62 0;			/* No iterative-direct algorithm */
	iparm[4] = 0;			/* No user fill-in reducing permutation */
	iparm[5] = 0;			/* Write solution into x */
	iparm[6] = 0;			/* Not in use */
	iparm[7] = 0;			/* Max numbers of iterative refinement steps */
	iparm[8] = 0;			/* Not in use */
	iparm[9] = 13;			/* Perturb the pivot elements with 1E-13 */
	iparm[10] = 1;			/* Use nonsymmetric permutation and scaling MPS */
	iparm[11] = 0;			/* Not in use */
	iparm[12] = 0;			/* Maximum weighted matching algorithm is switched-off (default for symmetric). Try iparm[12] = 1 in case of inappropriate accuracy */
	iparm[13] = 0;			/* Output: Number of perturbed pivots */
	iparm[14] = 0;			/* Not in use */
	iparm[15] = 0;			/* Not in use */
	iparm[16] = 0;			/* Not in use */
	iparm[17] = -1;			/* Output: Number of nonzeros in the factor LU */
	iparm[18] = -1;			/* Output: Mflops for LU factorization */
	iparm[19] = 0;			/* Output: Numbers of CG Iterations */

	iparm[59] = 0;			// OOC mode

	maxfct = 1;				/* Maximum number of numerical factorizations. */
	mnum = 1;				/* Which factorization to use. */
#ifdef DEBUG
	msglvl = 1;				/* Print statistical information in file */
#else
	msglvl = 0;
#endif
	error = 0;				/* Initialize error flag */
  msglvl=0;

	/* -------------------------------------------------------------------- */
	/* .. Initialize the internal solver memory pointer. This is only */
	/* necessary for the FIRST call of the PARDISO solver. */
	/* -------------------------------------------------------------------- */
	for (int i = 0; i < 64; i++)
	{
		pt[i] = 0;
	}

	m_nRhs		 = 1;
	m_factorized = 0;
}

SparseSolver::~SparseSolver() {

		this->Clear();

}

void SparseSolver::Clear() {

	if ( initialized == true )
	{
		double ddum;			/* Double dummy */
		MKL_INT idum;			/* Integer dummy. */
		MKL_INT nRhs = 1;

		/* -------------------------------------------------------------------- */
		/* .. Termination and release of memory. */
		/* -------------------------------------------------------------------- */
		phase = -1;			/* Release internal memory. */
		PARDISO (pt, &maxfct, &mnum, &mtype, &phase,
			&rows, &ddum, CSR_I_row_indices, CSR_J_col_indices, &idum, &nRhs,
			iparm, &msglvl, &ddum, &ddum, &error);

		initialized = false;

	}

	if (import_with_copy) {

		if (CSR_I_row_indices_size > 0)     delete [] CSR_I_row_indices;
		if (CSR_J_col_indices_size > 0)		delete [] CSR_J_col_indices;
		if (CSR_V_values_size > 0)			delete [] CSR_V_values;
	}

		CSR_I_row_indices_size = 0;
		CSR_J_col_indices_size = 0;
		CSR_V_values_size      = 0;

}

void SparseSolver::ImportMatrix(SparseMatrix & A) {

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

	import_with_copy = true;
}

void SparseSolver::ImportMatrix_wo_Copy(SparseMatrix & A) {

	rows	= A.rows;
	cols	= A.cols;
	nnz		= A.nnz;
	m_Kplus_size = A.rows;

	CSR_I_row_indices_size = A.CSR_I_row_indices.size();
	CSR_J_col_indices_size = A.CSR_J_col_indices.size();
	CSR_V_values_size	   = A.CSR_V_values.size();

	CSR_I_row_indices = &A.CSR_I_row_indices[0];
	CSR_J_col_indices = &A.CSR_J_col_indices[0];
	CSR_V_values	  = &A.CSR_V_values[0];

	import_with_copy = false;

}

void SparseSolver::SetThreaded() {

	/* Numbers of processors, value of OMP_NUM_THREADS */
	int num_procs;
	char * var = getenv("SOLVER_NUM_THREADS");
    if(var != NULL)
    	sscanf( var, "%d", &num_procs );
	else {
    	printf("Set environment SOLVER_NUM_THREADS to 1");
        exit(1);
	}

    iparm[2]  = num_procs;
}

void SparseSolver::Factorization() {

	double ddum;			/* Double dummy */
	MKL_INT idum;			/* Integer dummy. */

	/* -------------------------------------------------------------------- */
	/* .. Reordering and Symbolic Factorization. This step also allocates */
	/* all memory that is necessary for the factorization. */
	/* -------------------------------------------------------------------- */
	phase = 11;
	PARDISO (pt, &maxfct, &mnum, &mtype, &phase,
		&rows, CSR_V_values, CSR_I_row_indices, CSR_J_col_indices, &idum, &m_nRhs, iparm, &msglvl, &ddum, &ddum, &error);
	
	if (error != 0)
	{
		printf ("\nERROR during symbolic factorization: %d", error);
		exit (1);
	} else {
		initialized = true;
	}

#ifdef DEBUG
	printf ("\nReordering completed ... ");
	printf ("\nNumber of nonzeros in factors = %d", iparm[17]);
	printf ("\nNumber of factorization MFLOPS = %d", iparm[18]);
#endif

	/* -------------------------------------------------------------------- */
	/* .. Numerical factorization. */
	/* -------------------------------------------------------------------- */
	phase = 22;
	PARDISO (pt, &maxfct, &mnum, &mtype, &phase,
		&rows, CSR_V_values, CSR_I_row_indices, CSR_J_col_indices, &idum, &m_nRhs, iparm, &msglvl, &ddum, &ddum, &error);
	if (error != 0)
	{
		printf ("\nERROR during numerical factorization: %d", error);
		//exit (2);
	} else {
		m_factorized = 1;
	}

#ifdef DEBUG
	printf ("\nFactorization completed ... ");
#endif

	//TODO:
	tmp_sol.resize(m_Kplus_size); // - POZOR mozna se musi odkomentovat kvuli alokaci tmp_sol
}

void SparseSolver::Solve( SEQ_VECTOR <double> & rhs_sol) {

	if (!initialized)
		Factorization();

	double ddum   = 0;			/* Double dummy */
	MKL_INT idum  = 0;			/* Integer dummy. */
	MKL_INT n_rhs = 1;

	/* -------------------------------------------------------------------- */
	/* .. Back substitution and iterative refinement. */
	/* -------------------------------------------------------------------- */
	phase = 33;

	MKL_INT ip5backup = iparm[5];

	//iparm[24] = 1;		// Parallel forward/backward solve control. - 1 - Intel MKL PARDISO uses the sequential forward and backward solve.

	iparm[5] = 1;			// The solver stores the solution on the right-hand side b.
	PARDISO (pt, &maxfct, &mnum, &mtype, &phase,
		&rows, CSR_V_values, CSR_I_row_indices, CSR_J_col_indices, &idum, &n_rhs, iparm, &msglvl, &rhs_sol[0], &tmp_sol[0], &error);
	iparm[5] = ip5backup; 

	if (error != 0)
	{
		printf ("\nERROR during solution: %d", error);
		exit (3);
	}

	if (!keep_factors) {
		/* -------------------------------------------------------------------- */
		/* .. Termination and release of memory. */
		/* -------------------------------------------------------------------- */
		phase = -1;			/* Release internal memory. */
		MKL_INT nRhs = 1;
		PARDISO (pt, &maxfct, &mnum, &mtype, &phase,
				&rows, &ddum, CSR_I_row_indices, CSR_J_col_indices, &idum, &nRhs,
				iparm, &msglvl, &ddum, &ddum, &error);
		initialized = false;
		if (MPIrank == 0) printf(".");
	}

}

void SparseSolver::Solve( SEQ_VECTOR <double> & rhs, SEQ_VECTOR <double> & sol, MKL_INT n_rhs) {

	if (!initialized)
		Factorization();

	double ddum  = 0;			/* Double dummy */
	MKL_INT idum = 0;			/* Integer dummy. */

	/* -------------------------------------------------------------------- */
	/* .. Back substitution and iterative refinement. */
	/* -------------------------------------------------------------------- */
	phase = 33;
	//iparm[7] = 2;			/* Max numbers of iterative refinement steps. */

	PARDISO (pt, &maxfct, &mnum, &mtype, &phase,
		&rows, CSR_V_values, CSR_I_row_indices, CSR_J_col_indices, &idum, &n_rhs, iparm, &msglvl, &rhs[0], &sol[0], &error);
	if (error != 0)
	{
		printf ("\nERROR during solution: %d", error);
		exit (3);
	}

	if (!keep_factors) {
		/* -------------------------------------------------------------------- */
		/* .. Termination and release of memory. */
		/* -------------------------------------------------------------------- */
		phase = -1;			/* Release internal memory. */
		MKL_INT nRhs = 1;
		PARDISO (pt, &maxfct, &mnum, &mtype, &phase,
				&rows, &ddum, CSR_I_row_indices, CSR_J_col_indices, &idum, &nRhs,
				iparm, &msglvl, &ddum, &ddum, &error);
		initialized = false;
		if (MPIrank == 0) printf(".");
	}



}

void SparseSolver::Solve( SEQ_VECTOR <double> & rhs, SEQ_VECTOR <double> & sol, MKL_INT rhs_start_index, MKL_INT sol_start_index) {

	if (!initialized)
		Factorization();

	double ddum  = 0;			/* Double dummy */
	MKL_INT idum = 0;			/* Integer dummy. */

	/* -------------------------------------------------------------------- */
	/* .. Back substitution and iterative refinement. */
	/* -------------------------------------------------------------------- */
	phase = 33;
	//iparm[7] = 2;			/* Max numbers of iterative refinement steps. */

	PARDISO (pt, &maxfct, &mnum, &mtype, &phase,
		&rows, CSR_V_values, CSR_I_row_indices, CSR_J_col_indices, &idum, &m_nRhs, iparm, &msglvl, &rhs[rhs_start_index], &sol[sol_start_index], &error);
	if (error != 0)
	{
		printf ("\nERROR during solution: %d", error);
		exit (3);
	}


	if (!keep_factors) {
		/* -------------------------------------------------------------------- */
		/* .. Termination and release of memory. */
		/* -------------------------------------------------------------------- */
		phase = -1;			/* Release internal memory. */
		MKL_INT nRhs = 1;
		PARDISO (pt, &maxfct, &mnum, &mtype, &phase,
				&rows, &ddum, CSR_I_row_indices, CSR_J_col_indices, &idum, &nRhs,
				iparm, &msglvl, &ddum, &ddum, &error);
		initialized = false;
		if (MPIrank == 0) printf(".");
	}



}


void SparseSolver::SolveMat_Sparse( SparseMatrix & A) {
	SolveMat_Sparse(A, A);
};

void SparseSolver::SolveMat_Sparse( SparseMatrix & A_in, SparseMatrix & B_out) {
	SolveMat_Sparse(A_in, B_out, 'T');
};

void SparseSolver::SolveMat_Sparse( SparseMatrix & A_in, SparseMatrix & B_out, char T_for_input_matrix_is_transposed_N_input_matrix_is_NOT_transposed ) {

	if (!initialized)
		Factorization();

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
		phase = -1;			/* Release internal memory. */
		MKL_INT nRhs = 1;
		double ddum;			/* Double dummy */
		MKL_INT idum;			/* Integer dummy. */
		PARDISO (pt, &maxfct, &mnum, &mtype, &phase,
				&rows, &ddum, CSR_I_row_indices, CSR_J_col_indices, &idum, &nRhs,
				iparm, &msglvl, &ddum, &ddum, &error);
		initialized = false;
	}
}


void SparseSolver::SolveMat_Dense( SparseMatrix & A ) {
	SolveMat_Dense(A, A);
}

void SparseSolver::SolveMat_Dense( SparseMatrix & A_in, SparseMatrix & B_out ) {

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
void SparseSolver::SolveMatF( SparseMatrix & A_in, SparseMatrix & B_out, bool isThreaded ) {

	/* Internal solver memory pointer pt, */
	/* 32-bit: int pt[64]; 64-bit: long int pt[64] */
	/* or void *pt[64] should be OK on both architectures */
	void *pt[64];

	/* Pardiso control parameters. */
	MKL_INT iparm[64];
	MKL_INT maxfct, mnum, phase, error; //, msglvl;
	/* Auxiliary variables. */
	MKL_INT i;
	double ddum;			/* Double dummy */
	MKL_INT idum;			/* Integer dummy. */

	/* -------------------------------------------------------------------- */
	/* .. Setup Pardiso control parameters. */
	/* -------------------------------------------------------------------- */
	for (i = 0; i < 64; i++) {
		iparm[i] = 0;
	}

	MKL_INT mtype = 2;

	iparm[0] = 1;		/* No solver default */
	iparm[1] = 2;		/* Fill-in reordering from METIS */
						/* Numbers of processors, value of OMP_NUM_THREADS */
	//iparm[2] = 8;		/* Not used in MKL PARDISO */ // TODO: zjistit co to je pro MKL to bylo 0


	if (isThreaded) {
		/* Numbers of processors, value of OMP_NUM_THREADS */
		int num_procs;
		char * var = getenv("SOLVER_NUM_THREADS");
	    if(var != NULL)
	    	sscanf( var, "%d", &num_procs );
		else {
	    	printf("Set environment SOLVER_NUM_THREADS to 1");
	        exit(1);
		}

	    iparm[2] = num_procs;
	} else {
		iparm[2] = 1;
	}



	iparm[3] = 0;		/* No iterative-direct algorithm */
	iparm[4] = 0;		/* No user fill-in reducing permutation */
	iparm[5] = 0;		/* Write solution into x */
	iparm[6] = 0;		/* Not in use */
	iparm[7] = 0;		/* Max numbers of iterative refinement steps */
	iparm[8] = 0;		/* Not in use */
	iparm[9] = 13;		/* Perturb the pivot elements with 1E-13 */
	iparm[10] = 0;		/* Use nonsymmetric permutation and scaling MPS */
	iparm[11] = 0;		/* Not in use */
	iparm[12] = 0;		/* Maximum weighted matching algorithm is switched-off */
						/* (default for symmetric). Try iparm[12] = 1 in case of inappropriate accuracy */
	iparm[13] = 0;		/* Output: Number of perturbed pivots */
	iparm[14] = 0;		/* Not in use */
	iparm[15] = 0;		/* Not in use */
	iparm[16] = 0;		/* Not in use */
	iparm[17] = -1;		/* Output: Number of nonzeros in the factor LU */
	iparm[18] = -1;		/* Output: Mflops for LU factorization */
	iparm[19] = 0;		/* Output: Numbers of CG Iterations */

	maxfct = 1;			/* Maximum number of numerical factorizations. */
	mnum   = 1;			/* Which factorization to use. */
	//msglvl = 0;			/* Supress printing statistical information */
	error  = 0;			/* Initialize error flag */

	/* -------------------------------------------------------------------- */
	/* .. Initialize the internal solver memory pointer. This is only */
	/* necessary for the FIRST call of the PARDISO solver. */
	/* -------------------------------------------------------------------- */

	for (i = 0; i < 64; i++) pt[i] = 0;

	MKL_INT job[8];

	MKL_INT m		= A_in.rows;
	MKL_INT n		= A_in.cols;
	MKL_INT nRhs	= A_in.cols;
	MKL_INT lda     = m;
	MKL_INT info;

	SEQ_VECTOR<double>  sol  (m * n, 0);

	bool clear_dense = false;

	//TODO: na konci se musi smazat !!
	if (A_in.dense_values.size() == 0) {
		A_in.ConvertCSRToDense(0);
		clear_dense = true;
	}

	SEQ_VECTOR<MKL_INT> perm (A_in.dense_values.size() , 0);
	for (MKL_INT ii = 0; ii < A_in.dense_values.size(); ii++)
		if (A_in.dense_values[ii] != 0.0)
			perm[ii] = 1;

	phase = 13;
	PARDISO (pt, &maxfct, &mnum, &mtype, &phase,
		    &rows, CSR_V_values, CSR_I_row_indices, CSR_J_col_indices, &perm[0], &nRhs, iparm, &msglvl, &A_in.dense_values[0], &sol[0], &error);


	if (error != 0)
	{
		printf ("\nERROR during the solution of the system : %d", error);
		exit (1);
	} else {
		initialized = true;
	}

	MKL_INT x = 0;

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

	B_out.CSR_I_row_indices.resize(m + 1);
	B_out.CSR_J_col_indices.resize(1);
	B_out.CSR_V_values.     resize(1);

	mkl_ddnscsr (
		job,
		&m, &n,
		&sol[0], &lda,
		&B_out.CSR_V_values[0], &B_out.CSR_J_col_indices[0], &B_out.CSR_I_row_indices[0],
		&info);

	// Convert solution matrix (SOL) to sparse format - convert step
	MKL_INT nnzmax = B_out.CSR_I_row_indices[m];//-1; POZOR

	B_out.CSR_J_col_indices.resize(nnzmax);
	B_out.CSR_V_values.     resize(nnzmax);

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

	SEQ_VECTOR<double>().swap( sol );

	/* -------------------------------------------------------------------- */
	/* .. Termination and release of memory. */
	/* -------------------------------------------------------------------- */
	phase = -1;			/* Release internal memory. */
	PARDISO (pt, &maxfct, &mnum, &mtype, &phase,
		&rows, &ddum, CSR_I_row_indices, CSR_J_col_indices, &idum, &nRhs,
		iparm, &msglvl, &ddum, &ddum, &error);

	//SEQ_VECTOR<double>().swap( A_in.dense_values );

	initialized = false;

	if (clear_dense) {
		SEQ_VECTOR<double>().swap( A_in.dense_values );
	}
}

void SparseSolver::Create_SC( SparseMatrix & SC_out, MKL_INT sc_size, bool isThreaded ) {

	/* Internal solver memory pointer pt, */
	/* 32-bit: int pt[64]; 64-bit: long int pt[64] */
	/* or void *pt[64] should be OK on both architectures */
	void *pt[64];

	/* Pardiso control parameters. */
	MKL_INT 	iparm[64];
	double  dparm[65];
	MKL_INT 	maxfct, mnum, phase, error;

	/* Auxiliary variables. */
	MKL_INT 	i;
	double 	ddum;			/* Double dummy */
	MKL_INT 	idum;			/* Integer dummy. */
	MKL_INT 	solver;

	/* -------------------------------------------------------------------- */
	/* .. Setup Pardiso control parameters. */
	/* -------------------------------------------------------------------- */
	for (i = 0; i < 64; i++) {
		iparm[i] = 0;
	}

	/* -------------------------------------------------------------------- */
	/* .. Initialize the internal solver memory pointer. This is only */
	/* necessary for the FIRST call of the PARDISO solver. */
	/* -------------------------------------------------------------------- */
	for (i = 0; i < 64; i++)
		pt[i] = 0;

	MKL_INT 	mtype = 2;

	/* Numbers of processors, value of OMP_NUM_THREADS */
	if (isThreaded) {
		/* Numbers of processors, value of OMP_NUM_THREADS */
		int num_procs;
		char * var = getenv("SOLVER_NUM_THREADS");
	    if(var != NULL)
	    	sscanf( var, "%d", &num_procs );
		else {
	    	printf("Set environment SOLVER_NUM_THREADS to 1");
	        exit(1);
		}

	    iparm[2] = num_procs;
	} else {
		iparm[2] = 1;
	}

//	iparm[0] = 1;		/* No solver default */
//	iparm[1] = 2;		/* Fill-in reordering from METIS */
//	iparm[2]			/* Numbers of processors, value of OMP_NUM_THREADS */
//	iparm[2] = 8;		/* Not used in MKL PARDISO */
//	iparm[3] = 0;		/* No iterative-direct algorithm */
//	iparm[4] = 0;		/* No user fill-in reducing permutation */
//	iparm[5] = 0;		/* Write solution into x */
//	iparm[6] = 0;		/* Not in use */
//	iparm[7] = 0;		/* Max numbers of iterative refinement steps */
//	iparm[8] = 0;		/* Not in use */
//	iparm[9] = 13;		/* Perturb the pivot elements with 1E-13 */
//	iparm[10] = 0;		/* Use nonsymmetric permutation and scaling MPS */
//	iparm[11] = 0;		/* Not in use */
//	iparm[12] = 0;		/* Maximum weighted matching algorithm is switched-off */
//						/* (default for symmetric). Try iparm[12] = 1 in case of inappropriate accuracy */
//	iparm[13] = 0;		/* Output: Number of perturbed pivots */
//	iparm[14] = 0;		/* Not in use */
//	iparm[15] = 0;		/* Not in use */
//	iparm[16] = 0;		/* Not in use */
//	iparm[17] = -1;		/* Output: Number of nonzeros in the factor LU */
//	iparm[18] = -1;		/* Output: Mflops for LU factorization */
//	iparm[19] = 0;		/* Output: Numbers of CG Iterations */
//
//	maxfct = 1;			/* Maximum number of numerical factorizations. */
//	mnum   = 1;			/* Which factorization to use. */
//	//msglvl = 0;			/* Supress printing statistical information */
//	error  = 0;			/* Initialize error flag */




    iparm[1-1] = 1;         /* No solver default */
    iparm[2-1] = 2;         /* Fill-in reordering from METIS */
    iparm[10-1] = 8; //13   /* Perturb the pivot elements with 1E-13 */
    iparm[11-1] = 0;        /* Use nonsymmetric permutation and scaling MPS */
    iparm[13-1] = 0;        /* Maximum weighted matching algorithm is switched-off (default for symmetric). Try iparm[12] = 1 in case of inappropriate accuracy */
    iparm[14-1] = 0;        /* Output: Number of perturbed pivots */
    iparm[18-1] = -1;       /* Output: Number of nonzeros in the factor LU */
    iparm[19-1] = -1;       /* Output: Mflops for LU factorization */
    iparm[36-1] = 1;        /* Use Schur complement */

    maxfct = 1;           /* Maximum number of numerical factorizations. */
    mnum = 1;             /* Which factorization to use. */
    //msglvl = 1;           /* Print statistical information in file */
    error = 0;            /* Initialize error flag */

    /* -------------------------------------------------------------------- */
    /* .. Reordering and Symbolic Factorization. This step also allocates   */
    /* all memory that is necessary for the factorization.                  */
    /* -------------------------------------------------------------------- */

    std::vector <MKL_INT> perm (rows ,0);
    for (MKL_INT i = rows - sc_size; i < rows; i++)
    	perm[i] = 1;

    MKL_INT nrhs = 0;

//    phase = 11;
//    PARDISO (pt, &maxfct, &mnum, &mtype, &phase,
//        		&K_sc1.rows,
//				&K_sc1.CSR_V_values[0], &K_sc1.CSR_I_row_indices[0], &K_sc1.CSR_J_col_indices[0],
//				&perm[0], &nrhs, iparm, &msglvl, &ddum, &ddum, &error);
//
//    if ( error != 0 )
//    {
//    	printf ("\nERROR during symbolic factorization: %d", error);
//    	exit(1);
//    }


    /* -------------------------------------------------------------------- */
    /* .. Numerical factorization. */
    /* -------------------------------------------------------------------- */

	SC_out.dense_values.resize(sc_size * sc_size);

    phase = 12;
    PARDISO (pt, &maxfct, &mnum, &mtype, &phase,
			&rows,
			&CSR_V_values[0], &CSR_I_row_indices[0], &CSR_J_col_indices[0],
			&perm[0], &nrhs,
			iparm, &msglvl, &ddum, &SC_out.dense_values[0], &error);

    //for (MKL_INT i = 0; i < SC_out.dense_values.size(); i++)
    //	SC_out.dense_values[i] = (-1.0)*SC_out.dense_values[i];

    if ( error != 0 )
	{
		printf ("\nERROR during numerical factorization: %d", error);
		exit (2);
	} else {
		initialized = true;
	}

	/* -------------------------------------------------------------------- */
	/* .. Termination and release of memory. */
	/* -------------------------------------------------------------------- */
	phase = -1;           /* Release internal memory. */
	PARDISO (pt, &maxfct, &mnum, &mtype, &phase,
			&rows, &ddum, &CSR_I_row_indices[0], &CSR_J_col_indices[0], &idum, &nrhs,
			 iparm, &msglvl, &ddum, &ddum, &error);

	initialized = false;

    /* -------------------------------------------------------------------- */
    /* ..  allocate memory for the Schur-complement and copy it there.      */
    /* -------------------------------------------------------------------- */
    MKL_INT nonzeros_S = iparm[38];

    SC_out.cols = sc_size;
    SC_out.rows = sc_size;
    SC_out.type = 'G';

    SC_out.ConvertDenseToCSR(1);
    SC_out.RemoveLower();
    SC_out.type = 'S';

    if (msglvl == 1)
    	SpyText(SC_out);

//    if (generate_symmetric_sc_1_generate_general_sc_0 == 1) {
//    	SC_out.RemoveLower();
//    }

}

void SparseSolver::Create_SC_w_Mat( SparseMatrix & K_in, SparseMatrix & B_in, SparseMatrix & SC_out,
								    bool isThreaded, MKL_INT generate_symmetric_sc_1_generate_general_sc_0 ) {

	//int msglvl = 0;

	// *** Prepare matrix

	SparseMatrix K_sc1;
	SparseMatrix Sc_eye;
	SparseMatrix K_b_tmp;

	K_b_tmp = B_in;
	K_b_tmp.MatTranspose();

	Sc_eye.CreateEye(K_b_tmp.rows, 0.0, 0, K_b_tmp.cols);

	K_sc1 = K_in;
	K_sc1.MatTranspose();
	K_sc1.MatAppend(K_b_tmp);
	K_sc1.MatTranspose();
	K_sc1.MatAppend(Sc_eye);


	// *** END - Prepare matrix




	/* Internal solver memory pointer pt, */
	/* 32-bit: int pt[64]; 64-bit: long int pt[64] */
	/* or void *pt[64] should be OK on both architectures */
	void *pt[64];

	/* Pardiso control parameters. */
	MKL_INT 	iparm[64];
	double  	dparm[65];
	MKL_INT 	maxfct, mnum, phase, error;

	/* Auxiliary variables. */
	MKL_INT 	i;
	double 		ddum;			/* Double dummy */
	MKL_INT 	idum;			/* Integer dummy. */
	MKL_INT 	solver;

	/* -------------------------------------------------------------------- */
	/* .. Setup Pardiso control parameters. */
	/* -------------------------------------------------------------------- */
	for (i = 0; i < 64; i++) {
		iparm[i] = 0;
	}

	/* -------------------------------------------------------------------- */
	/* .. Initialize the internal solver memory pointer. This is only */
	/* necessary for the FIRST call of the PARDISO solver. */
	/* -------------------------------------------------------------------- */
	for (i = 0; i < 64; i++)
		pt[i] = 0;

	MKL_INT 	mtype = 2;

	/* Numbers of processors, value of OMP_NUM_THREADS */
	if (isThreaded) {
		/* Numbers of processors, value of OMP_NUM_THREADS */
		MKL_INT num_procs;
		char * var = getenv("SOLVER_NUM_THREADS");
	    if(var != NULL)
	    	sscanf( var, "%d", &num_procs );
		else {
	    	printf("Set environment SOLVER_NUM_THREADS to 1");
	        exit(1);
		}

	    iparm[2] = num_procs;
	} else {
		iparm[2] = 1;
	}

//	iparm[0] = 1;		/* No solver default */
//	iparm[1] = 2;		/* Fill-in reordering from METIS */
//	iparm[2]			/* Numbers of processors, value of OMP_NUM_THREADS */
//	iparm[2] = 8;		/* Not used in MKL PARDISO */
//	iparm[3] = 0;		/* No iterative-direct algorithm */
//	iparm[4] = 0;		/* No user fill-in reducing permutation */
//	iparm[5] = 0;		/* Write solution into x */
//	iparm[6] = 0;		/* Not in use */
//	iparm[7] = 0;		/* Max numbers of iterative refinement steps */
//	iparm[8] = 0;		/* Not in use */
//	iparm[9] = 13;		/* Perturb the pivot elements with 1E-13 */
//	iparm[10] = 0;		/* Use nonsymmetric permutation and scaling MPS */
//	iparm[11] = 0;		/* Not in use */
//	iparm[12] = 0;		/* Maximum weighted matching algorithm is switched-off */
//						/* (default for symmetric). Try iparm[12] = 1 in case of inappropriate accuracy */
//	iparm[13] = 0;		/* Output: Number of perturbed pivots */
//	iparm[14] = 0;		/* Not in use */
//	iparm[15] = 0;		/* Not in use */
//	iparm[16] = 0;		/* Not in use */
//	iparm[17] = -1;		/* Output: Number of nonzeros in the factor LU */
//	iparm[18] = -1;		/* Output: Mflops for LU factorization */
//	iparm[19] = 0;		/* Output: Numbers of CG Iterations */
//
//	maxfct = 1;			/* Maximum number of numerical factorizations. */
//	mnum   = 1;			/* Which factorization to use. */
//	//msglvl = 0;			/* Supress printing statistical information */
//	error  = 0;			/* Initialize error flag */




    iparm[1-1] = 1;         /* No solver default */
    iparm[2-1] = 2;         /* Fill-in reordering from METIS */
    iparm[10-1] = 8; //13   /* Perturb the pivot elements with 1E-13 */
    iparm[11-1] = 0;        /* Use nonsymmetric permutation and scaling MPS */
    iparm[13-1] = 0;        /* Maximum weighted matching algorithm is switched-off (default for symmetric). Try iparm[12] = 1 in case of inappropriate accuracy */
    iparm[14-1] = 0;        /* Output: Number of perturbed pivots */
    iparm[18-1] = -1;       /* Output: Number of nonzeros in the factor LU */
    iparm[19-1] = -1;       /* Output: Mflops for LU factorization */
    iparm[36-1] = 1;        /* Use Schur complement */

    maxfct = 1;           /* Maximum number of numerical factorizations. */
    mnum = 1;             /* Which factorization to use. */
    //msglvl = 1;           /* Print statistical information in file */
    error = 0;            /* Initialize error flag */

    /* -------------------------------------------------------------------- */
    /* .. Reordering and Symbolic Factorization. This step also allocates   */
    /* all memory that is necessary for the factorization.                  */
    /* -------------------------------------------------------------------- */

    std::vector <MKL_INT> perm (K_sc1.rows,0);
    for (MKL_INT i = K_in.rows; i < K_sc1.rows; i++)
    	perm[i] = 1;

    MKL_INT nrhs = 0;

//    phase = 11;
//    PARDISO (pt, &maxfct, &mnum, &mtype, &phase,
//        		&K_sc1.rows,
//				&K_sc1.CSR_V_values[0], &K_sc1.CSR_I_row_indices[0], &K_sc1.CSR_J_col_indices[0],
//				&perm[0], &nrhs, iparm, &msglvl, &ddum, &ddum, &error);
//
//    if ( error != 0 )
//    {
//    	printf ("\nERROR during symbolic factorization: %d", error);
//    	exit(1);
//    }


    /* -------------------------------------------------------------------- */
    /* .. Numerical factorization. */
    /* -------------------------------------------------------------------- */

	SC_out.dense_values.resize(K_b_tmp.rows * K_b_tmp.rows);

    phase = 12;
    PARDISO (pt, &maxfct, &mnum, &mtype, &phase,
			&K_sc1.rows,
			&K_sc1.CSR_V_values[0], &K_sc1.CSR_I_row_indices[0], &K_sc1.CSR_J_col_indices[0],
			&perm[0], &nrhs,
			iparm, &msglvl, &ddum, &SC_out.dense_values[0], &error);

    for (MKL_INT i = 0; i < SC_out.dense_values.size(); i++)
    	SC_out.dense_values[i] = (-1.0)*SC_out.dense_values[i];

    if ( error != 0 )
	{
		printf ("\nERROR during numerical factorization: %d", error);
		exit (2);
	} else {
		initialized = true;
	}

	/* -------------------------------------------------------------------- */
	/* .. Termination and release of memory. */
	/* -------------------------------------------------------------------- */
	phase = -1;           /* Release internal memory. */
	PARDISO (pt, &maxfct, &mnum, &mtype, &phase,
			&K_sc1.rows, &ddum, &K_sc1.CSR_I_row_indices[0], &K_sc1.CSR_J_col_indices[0], &idum, &nrhs,
			 iparm, &msglvl, &ddum, &ddum, &error);

	initialized = false;

    /* -------------------------------------------------------------------- */
    /* ..  allocate memory for the Schur-complement and copy it there.      */
    /* -------------------------------------------------------------------- */
    MKL_INT nonzeros_S = iparm[38];

    SC_out.cols = K_b_tmp.rows;
    SC_out.rows = K_b_tmp.rows;
    SC_out.type = 'G';

//    SC_out.ConvertDenseToCSR(1);

//    if (msglvl == 1)
//    	SpyText(SC_out);

    if (generate_symmetric_sc_1_generate_general_sc_0 == 1) {
    	//SC_out.RemoveLower();
    	SC_out.RemoveLowerDense();
    	SC_out.type = 'S';
    }


//    // Finalize shape of the SC
//    if (generate_symmetric_sc_1_generate_general_sc_0 == 0) {
//
//		SC_out.type = 'G';
//
//		SparseMatrix SC_tmp;
//		SC_tmp = SC_out;
//		SC_tmp.SetDiagonalOfSymmetricMatrix(0.0);
//		SC_tmp.MatTranspose();
//
//		SC_out.MatAddInPlace(SC_tmp,'N',1.0);
//
//    }
//
//	SC_out.MatScale(-1.0);
//	//SC_out.ConvertCSRToDense(0);

}

void SparseSolver::Create_non_sym_SC_w_Mat( SparseMatrix & K_in, SparseMatrix & B1_in, SparseMatrix & B0_in, SparseMatrix & SC_out, bool isThreaded, MKL_INT generate_symmetric_sc_1_generate_general_sc_0 ) {

	//int msglvl = 0;

	//SpyText(K_in);
	//SpyText(B1_in);
	//SpyText(B0_in);


	SparseMatrix K;

	// Create "non-symmetric" K matrix;
	if (K_in.type = 'S') {
		K = K_in;
		K.SetDiagonalOfSymmetricMatrix(0.0);
		K.MatTranspose();
		K.MatAddInPlace(K_in, 'N', 1.0);
		K.type = 'G';
	} else {
		// Pozor - not optimal, zbytecna kopie matice K
		K = K_in;
	}


	//SpyText(K);

	// *** Prepare matrix

	SparseMatrix K_sc1;
	SparseMatrix Sc_eye;
	SparseMatrix K_b_tmp;

	K_b_tmp = B1_in;
	K_b_tmp.MatTranspose();

	Sc_eye.CreateEye(K_b_tmp.rows, 0.0, 0, K_b_tmp.cols);

	K.MatTranspose();
	K.MatAppend(K_b_tmp);
	K.MatTranspose();

    //SpyText(K);

	SparseMatrix K_b0_tmp;
	K_b0_tmp = B0_in;
	K_b0_tmp.MatTranspose();
	K_b0_tmp.ConvertToCOO(1);
	K_b0_tmp.rows = B1_in.cols;
	K_b0_tmp.cols = B1_in.cols + K_in.cols;
	K_b0_tmp.ConvertToCSRwithSort(1);

    //SpyText(K_b0_tmp);

	K_b0_tmp.MatAddInPlace(Sc_eye,'N',1.0);

    //SpyText(K_b0_tmp);

	K.MatAppend(K_b0_tmp);

    K_sc1 = K;

	// *** END - Prepare matrix
	//SpyText(K_sc1);


	/* Internal solver memory pointer pt, */
	/* 32-bit: int pt[64]; 64-bit: long int pt[64] */
	/* or void *pt[64] should be OK on both architectures */
	void *pt[64];

	/* Pardiso control parameters. */
	MKL_INT 	iparm[64];
	double  dparm[65];
	MKL_INT 	maxfct, mnum, phase, error;

	/* Auxiliary variables. */
	MKL_INT 	i;
	double 	ddum;			/* Double dummy */
	MKL_INT 	idum;			/* Integer dummy. */
	MKL_INT 	solver;

	/* -------------------------------------------------------------------- */
	/* .. Setup Pardiso control parameters. */
	/* -------------------------------------------------------------------- */
	for (i = 0; i < 64; i++) {
		iparm[i] = 0;
	}

	/* -------------------------------------------------------------------- */
	/* .. Initialize the internal solver memory pointer. This is only */
	/* necessary for the FIRST call of the PARDISO solver. */
	/* -------------------------------------------------------------------- */
	for (i = 0; i < 64; i++)
		pt[i] = 0;

	MKL_INT 	mtype = 11;

	/* Numbers of processors, value of OMP_NUM_THREADS */
	if (isThreaded) {
		/* Numbers of processors, value of OMP_NUM_THREADS */
		int num_procs;
		char * var = getenv("SOLVER_NUM_THREADS");
	    if(var != NULL)
	    	sscanf( var, "%d", &num_procs );
		else {
	    	printf("Set environment SOLVER_NUM_THREADS to 1");
	        exit(1);
		}

	    iparm[2] = num_procs;
	} else {
		iparm[2] = 1;
	}

//	iparm[0] = 1;		/* No solver default */
//	iparm[1] = 2;		/* Fill-in reordering from METIS */
//	iparm[2]			/* Numbers of processors, value of OMP_NUM_THREADS */
//	iparm[2] = 8;		/* Not used in MKL PARDISO */
//	iparm[3] = 0;		/* No iterative-direct algorithm */
//	iparm[4] = 0;		/* No user fill-in reducing permutation */
//	iparm[5] = 0;		/* Write solution into x */
//	iparm[6] = 0;		/* Not in use */
//	iparm[7] = 0;		/* Max numbers of iterative refinement steps */
//	iparm[8] = 0;		/* Not in use */
//	iparm[9] = 13;		/* Perturb the pivot elements with 1E-13 */
//	iparm[10] = 0;		/* Use nonsymmetric permutation and scaling MPS */
//	iparm[11] = 0;		/* Not in use */
//	iparm[12] = 0;		/* Maximum weighted matching algorithm is switched-off */
//						/* (default for symmetric). Try iparm[12] = 1 in case of inappropriate accuracy */
//	iparm[13] = 0;		/* Output: Number of perturbed pivots */
//	iparm[14] = 0;		/* Not in use */
//	iparm[15] = 0;		/* Not in use */
//	iparm[16] = 0;		/* Not in use */
//	iparm[17] = -1;		/* Output: Number of nonzeros in the factor LU */
//	iparm[18] = -1;		/* Output: Mflops for LU factorization */
//	iparm[19] = 0;		/* Output: Numbers of CG Iterations */
//
//	maxfct = 1;			/* Maximum number of numerical factorizations. */
//	mnum   = 1;			/* Which factorization to use. */
//	//msglvl = 0;			/* Supress printing statistical information */
//	error  = 0;			/* Initialize error flag */




    iparm[1-1] = 1;         /* No solver default */
    iparm[2-1] = 2;         /* Fill-in reordering from METIS */
    iparm[10-1] = 8; //13   /* Perturb the pivot elements with 1E-13 */
    iparm[11-1] = 0;        /* Use nonsymmetric permutation and scaling MPS */
    iparm[13-1] = 0;        /* Maximum weighted matching algorithm is switched-off (default for symmetric). Try iparm[12] = 1 in case of inappropriate accuracy */
    iparm[14-1] = 0;        /* Output: Number of perturbed pivots */
    iparm[18-1] = -1;       /* Output: Number of nonzeros in the factor LU */
    iparm[19-1] = -1;       /* Output: Mflops for LU factorization */
    iparm[36-1] = 1;        /* Use Schur complement */

    maxfct = 1;           /* Maximum number of numerical factorizations. */
    mnum = 1;             /* Which factorization to use. */
    //msglvl = 1;           /* Print statistical information in file */
    error = 0;            /* Initialize error flag */

    /* -------------------------------------------------------------------- */
    /* .. Reordering and Symbolic Factorization. This step also allocates   */
    /* all memory that is necessary for the factorization.                  */
    /* -------------------------------------------------------------------- */

    std::vector <MKL_INT> perm (K_sc1.rows,0);
    for (MKL_INT i = K_in.rows; i < K_sc1.rows; i++)
    	perm[i] = 1;

    MKL_INT nrhs = 0;

//    phase = 11;
//    PARDISO (pt, &maxfct, &mnum, &mtype, &phase,
//        		&K_sc1.rows,
//				&K_sc1.CSR_V_values[0], &K_sc1.CSR_I_row_indices[0], &K_sc1.CSR_J_col_indices[0],
//				&perm[0], &nrhs, iparm, &msglvl, &ddum, &ddum, &error);
//
//    if ( error != 0 )
//    {
//    	printf ("\nERROR during symbolic factorization: %d", error);
//    	exit(1);
//    }


    /* -------------------------------------------------------------------- */
    /* .. Numerical factorization. */
    /* -------------------------------------------------------------------- */

	SC_out.dense_values.resize(K_b_tmp.rows * K_b_tmp.rows);

    phase = 12;
    PARDISO (pt, &maxfct, &mnum, &mtype, &phase,
			&K_sc1.rows,
			&K_sc1.CSR_V_values[0], &K_sc1.CSR_I_row_indices[0], &K_sc1.CSR_J_col_indices[0],
			&perm[0], &nrhs,
			iparm, &msglvl, &ddum, &SC_out.dense_values[0], &error);

    for (MKL_INT i = 0; i < SC_out.dense_values.size(); i++)
    	SC_out.dense_values[i] = (-1.0)*SC_out.dense_values[i];

    if ( error != 0 )
	{
		printf ("\nERROR during numerical factorization: %d", error);
		exit (2);
	} else {
		initialized = true;
	}

	/* -------------------------------------------------------------------- */
	/* .. Termination and release of memory. */
	/* -------------------------------------------------------------------- */
	phase = -1;           /* Release internal memory. */
	PARDISO (pt, &maxfct, &mnum, &mtype, &phase,
			&K_sc1.rows, &ddum, &K_sc1.CSR_I_row_indices[0], &K_sc1.CSR_J_col_indices[0], &idum, &nrhs,
			 iparm, &msglvl, &ddum, &ddum, &error);

	initialized = false;

    /* -------------------------------------------------------------------- */
    /* ..  allocate memory for the Schur-complement and copy it there.      */
    /* -------------------------------------------------------------------- */

    SC_out.cols = K_b_tmp.rows;
    SC_out.rows = B0_in.cols;
    SC_out.type = 'G';
    SC_out.dense_values.resize( SC_out.cols * SC_out.rows);

    //SC_out.ConvertDenseToCSR(1);

}


#include <stdio.h>
#include "mkl_rci.h"
#include "mkl_blas.h"
#include "mkl_spblas.h"
#include "mkl_service.h"

void SparseSolver::SolveCG(SparseMatrix & A_in, SEQ_VECTOR <double> & rhs_sol) {

	MKL_INT size = A_in.rows;
	SEQ_VECTOR <double> sol (size, 0);

	SolveCG(A_in, rhs_sol, sol);

	rhs_sol = sol;
}

void SparseSolver::SolveCG(SparseMatrix & A_in, SEQ_VECTOR <double> & rhs_in, SEQ_VECTOR <double> & sol) {




	  /*---------------------------------------------------------------------------   */
	  /* Define arrays for the upper triangle of the coefficient matrix and rhs vector */
	  /* Compressed sparse row storage is used for sparse representation              */
	  /*---------------------------------------------------------------------------   */
	  MKL_INT rci_request, itercount, i;

	  //MKL_INT expected_itercount = 8;
	  MKL_INT n = A_in.rows; //rhs_in.size();


	  /* Fill all arrays containing matrix data. */
	  MKL_INT * ia  = &A_in.CSR_I_row_indices[0];
	  MKL_INT * ja  = &A_in.CSR_J_col_indices[0];
	  double  * a   = &A_in.CSR_V_values[0];


//	  MKL_INT * ia  = CSR_I_row_indices;
//	  MKL_INT * ja  = CSR_J_col_indices;
//	  double  * a   = CSR_V_values;

	  /*---------------------------------------------------------------------------*/
	  /* Allocate storage for the solver ?par and temporary storage tmp            */
	  /*---------------------------------------------------------------------------*/
	  MKL_INT length = 128;
	  MKL_INT ipar[128];
	  double dpar[128];
	  char matdes[3];
	  double one = 1.E0;

	  SEQ_VECTOR <double> tmp_vec (4 * n, 0);
	  double * tmp;
	  tmp = &tmp_vec[0];

	  /*---------------------------------------------------------------------------*/
	  /* Some additional variables to use with the RCI (P)CG solver                */
	  /*---------------------------------------------------------------------------*/
	  double * solution = &sol[0];
      double * rhs = &rhs_in[0];


	  double * temp;
	  SEQ_VECTOR <double> temp_vec (n,0);
	  temp = &temp_vec[0];

	  double euclidean_norm;

	  char tr = 'u';
	  double eone = -1.E0;
	  MKL_INT ione = 1;

	  /*---------------------------------------------------------------------------*/
	  /* Initialize the initial guess                                              */
	  /*---------------------------------------------------------------------------*/
	  for (i = 0; i < n; i++)
	    {
	      solution[i] = 0.E0;
	    }

	  matdes[0] = 'd';
	  matdes[1] = 'l';
	  matdes[2] = 'n';
	  /*---------------------------------------------------------------------------*/
	  /* Initialize the solver                                                     */
	  /*---------------------------------------------------------------------------*/

	  dcg_init (&n, solution, rhs, &rci_request, ipar, dpar, tmp);
	  if (rci_request != 0)
	    goto failure;
	  /*---------------------------------------------------------------------------*/
	  /* Set the desired parameters:                                               */
	  /* LOGICAL parameters:                                                       */
	  /* -                                                                         */
	  /* INTEGER parameters:                                                       */
	  /* set the maximal number of iterations to 100                               */
	  /* DOUBLE parameters                                                         */
	  /* -                                                                         */
	  /*---------------------------------------------------------------------------*/
	  ipar[4] = 1000;
	  ipar[10]=1;
	  /*---------------------------------------------------------------------------*/
	  /* Check the correctness and consistency of the newly set parameters         */
	  /*---------------------------------------------------------------------------*/
	  dcg_check (&n, solution, rhs, &rci_request, ipar, dpar, tmp);
	  if (rci_request != 0)
	    goto failure;
	  /*---------------------------------------------------------------------------*/
	  /* Compute the solution by RCI (P)CG solver                                  */
	  /* Reverse Communications starts here                                        */
	  /*---------------------------------------------------------------------------*/
	rci:dcg (&n, solution, rhs, &rci_request, ipar, dpar, tmp);
	  /*---------------------------------------------------------------------------*/
	  /* If rci_request=0, then the solution was found according to the requested  */
	  /* stopping tests. In this case, this means that it was found after 100      */
	  /* iterations.                                                               */
	  /*---------------------------------------------------------------------------*/
	  if (rci_request == 0)
	    goto getsln;
	  /*---------------------------------------------------------------------------*/
	  /* If rci_request=1, then compute the vector A*TMP[0]                        */
	  /* and put the result in vector TMP[n]                                       */
	  /*---------------------------------------------------------------------------*/
	  if (rci_request == 1)
	    {
	      mkl_dcsrsymv (&tr, &n, a, ia, ja, tmp, &tmp[n]);
	      goto rci;
	    }
	  /*---------------------------------------------------------------------------*/
	  /* If rci_request=2, then do the user-defined stopping test: compute the     */
	  /* Euclidean norm of the actual residual using MKL routines and check if     */
	  /* it is less than 1.E-8                                                     */
	  /*---------------------------------------------------------------------------*/
	  if (rci_request == 2)
	    {
	      mkl_dcsrsymv (&tr, &n, a, ia, ja, solution, temp);
	      daxpy (&n, &eone, rhs, &ione, temp, &ione);
	      euclidean_norm = dnrm2 (&n, temp, &ione);
	      /*---------------------------------------------------------------------------*/
	      /* The solution has not been found yet according to the user-defined stopping */
	      /* test. Continue RCI (P)CG iterations.                                      */
	      /*---------------------------------------------------------------------------*/
	      if (euclidean_norm > 1.E-6)
	        goto rci;
	      /*---------------------------------------------------------------------------*/
	      /* The solution has been found according to the user-defined stopping test   */
	      /*---------------------------------------------------------------------------*/
	      else
	        goto getsln;
	    }
	  /*---------------------------------------------------------------------------*/
	  /* If rci_request=3, then compute apply the preconditioner matrix C_inverse  */
	  /* on vector tmp[2*n] and put the result in vector tmp[3*n]                  */
	  /*---------------------------------------------------------------------------*/
	  if (rci_request == 3)
	    {
	      mkl_dcsrsv (&matdes[2], &n, &one, matdes, a, ja, ia, &ia[1], &tmp[2 * n], &tmp[3 * n]);
	      goto rci;
	    }	  /*---------------------------------------------------------------------------*/
	  /* If rci_request=anything else, then dcg subroutine failed                  */
	  /* to compute the solution vector: solution[n]                               */
	  /*---------------------------------------------------------------------------*/
	  goto failure;
	  /*---------------------------------------------------------------------------*/
	  /* Reverse Communication ends here                                           */
	  /* Get the current iteration number into itercount                           */
	  /*---------------------------------------------------------------------------*/
	getsln:dcg_get (&n, solution, rhs, &rci_request, ipar, dpar, tmp, &itercount);

	std::cout << " " << itercount;

	  //printf ("\nNumber of iterations: %d\n", itercount);

	  /*-------------------------------------------------------------------------*/
	  /* Release internal MKL memory that might be used for computations         */
	  /* NOTE: It is important to call the routine below to avoid memory leaks   */
	  /* unless you disable MKL Memory Manager                                   */
	  /*-------------------------------------------------------------------------*/
	failure: ; //printf ("This example FAILED as the solver has returned the ERROR code %d", rci_request);
	  //MKL_Free_Buffers ();


}
