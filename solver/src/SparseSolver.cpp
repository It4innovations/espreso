#include "SparseSolver.h"

/* PARDISO prototype. */
extern "C" void pardisoinit (void   *, int    *,   int *, int *, double *, int *);
extern "C" void pardiso     (void   *, int    *,   int *, int *,    int *, int *,
                  double *, int    *,    int *, int *,   int *, int *,
                     int *, double *, double *, int *, double *);
extern "C" void pardiso_chkmatrix  (int *, int *, double *, int *, int *, int *);
extern "C" void pardiso_chkvec     (int *, int *, double *, int *);
extern "C" void pardiso_printstats (int *, int *, double *, int *, int *, int *,
                           double *, int *);
extern "C" void pardiso_get_schur  (void*, int*, int*, int*, double*, int*, int*);


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

	iparm[3] = 0;			/* No iterative-direct algorithm */
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

	iparm[59] = 1;			// OOC mode

	maxfct = 1;				/* Maximum number of numerical factorizations. */
	mnum = 1;				/* Which factorization to use. */
#ifdef DEBUG
	msglvl = 1;				/* Print statistical information in file */
#else
	msglvl = 0;
#endif
	error = 0;				/* Initialize error flag */

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
		int idum;			/* Integer dummy. */
		int nRhs = 1;

		/* -------------------------------------------------------------------- */
		/* .. Termination and release of memory. */
		/* -------------------------------------------------------------------- */
		phase = -1;			/* Release internal memory. */

		pardiso (pt, &maxfct, &mnum, &mtype, &phase,
			&rows, &ddum, CSR_I_row_indices, CSR_J_col_indices, &idum, &nRhs,
			iparm, &msglvl, &ddum, &ddum, &error, dparm);

		initialized = false;

	}

		if (CSR_I_row_indices_size > 0)     delete [] CSR_I_row_indices;
		if (CSR_J_col_indices_size > 0)		delete [] CSR_J_col_indices;
		if (CSR_V_values_size > 0)			delete [] CSR_V_values;

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

	CSR_I_row_indices = new int[CSR_I_row_indices_size];
	CSR_J_col_indices = new int[CSR_J_col_indices_size];
	CSR_V_values	  = new double  [CSR_V_values_size];

	if (sizeof(A.CSR_I_row_indices[0]) == sizeof(int)) {
		copy(A.CSR_I_row_indices.begin(), A.CSR_I_row_indices.end(), CSR_I_row_indices);
		copy(A.CSR_J_col_indices.begin(), A.CSR_J_col_indices.end(), CSR_J_col_indices);
		copy(A.CSR_V_values     .begin(), A.CSR_V_values     .end(), CSR_V_values);
	} else {
		for (int i = 0; i < A.CSR_I_row_indices.size(); i++)
			CSR_I_row_indices[i] = A.CSR_I_row_indices[i];

		for (int i = 0; i < A.CSR_J_col_indices.size(); i++)
			CSR_J_col_indices[i] = A.CSR_J_col_indices[i];

		for (int i = 0; i < A.CSR_V_values.size(); i++)
			CSR_V_values[i] = A.CSR_V_values[i];
	}
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
	int idum;				/* Integer dummy. */

	/* -------------------------------------------------------------------- */
	/* .. Reordering and Symbolic Factorization. This step also allocates */
	/* all memory that is necessary for the factorization. */
	/* -------------------------------------------------------------------- */
	phase = 11;
	pardiso (pt, &maxfct, &mnum, &mtype, &phase,
		&rows, CSR_V_values, CSR_I_row_indices, CSR_J_col_indices, &idum, &m_nRhs, iparm, &msglvl, &ddum, &ddum, &error, dparm);

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
	pardiso (pt, &maxfct, &mnum, &mtype, &phase,
		&rows, CSR_V_values, CSR_I_row_indices, CSR_J_col_indices, &idum, &m_nRhs, iparm, &msglvl, &ddum, &ddum, &error, dparm);
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

	double ddum  = 0;			/* Double dummy */
	int idum     = 0;			/* Integer dummy. */
	int n_rhs    = 1;

	/* -------------------------------------------------------------------- */
	/* .. Back substitution and iterative refinement. */
	/* -------------------------------------------------------------------- */
	phase = 33;

	int ip5backup = iparm[5];

	//iparm[24] = 1;		// Parallel forward/backward solve control. - 1 - Intel MKL PARDISO uses the sequential forward and backward solve.

	iparm[5] = 1;			// The solver stores the solution on the right-hand side b.
	pardiso (pt, &maxfct, &mnum, &mtype, &phase,
		&rows, CSR_V_values, CSR_I_row_indices, CSR_J_col_indices, &idum, &n_rhs, iparm, &msglvl, &rhs_sol[0], &tmp_sol[0], &error, dparm);
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
		int nRhs = 1;
		pardiso (pt, &maxfct, &mnum, &mtype, &phase,
			&rows, &ddum, CSR_I_row_indices, CSR_J_col_indices, &idum, &nRhs,
			iparm, &msglvl, &ddum, &ddum, &error, dparm);
		initialized = false;
		if (MPIrank == 0) printf(".");
	}

}

void SparseSolver::Solve( SEQ_VECTOR <double> & rhs, SEQ_VECTOR <double> & sol, int n_rhs) {

	if (!initialized)
		Factorization();

	double ddum  = 0;			/* Double dummy */
	int idum = 0;			/* Integer dummy. */

	/* -------------------------------------------------------------------- */
	/* .. Back substitution and iterative refinement. */
	/* -------------------------------------------------------------------- */
	phase = 33;
	//iparm[7] = 2;			/* Max numbers of iterative refinement steps. */

	pardiso (pt, &maxfct, &mnum, &mtype, &phase,
		&rows, CSR_V_values, CSR_I_row_indices, CSR_J_col_indices, &idum, &n_rhs, iparm, &msglvl, &rhs[0], &sol[0], &error, dparm);

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
		int nRhs = 1;
		pardiso (pt, &maxfct, &mnum, &mtype, &phase,
			&rows, &ddum, CSR_I_row_indices, CSR_J_col_indices, &idum, &nRhs,
			iparm, &msglvl, &ddum, &ddum, &error, dparm);
		initialized = false;
		if (MPIrank == 0) printf(".");
	}

}

void SparseSolver::Solve( SEQ_VECTOR <double> & rhs, SEQ_VECTOR <double> & sol, int rhs_start_index, int sol_start_index) {

	if (!initialized)
		Factorization();

	double ddum  = 0;			/* Double dummy */
	int idum = 0;			/* Integer dummy. */

	/* -------------------------------------------------------------------- */
	/* .. Back substitution and iterative refinement. */
	/* -------------------------------------------------------------------- */
	phase = 33;
	//iparm[7] = 2;			/* Max numbers of iterative refinement steps. */

	pardiso (pt, &maxfct, &mnum, &mtype, &phase,
		&rows, CSR_V_values, CSR_I_row_indices, CSR_J_col_indices, &idum, &m_nRhs, iparm, &msglvl, &rhs[rhs_start_index], &sol[sol_start_index], &error, dparm);
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
		int nRhs = 1;
		pardiso (pt, &maxfct, &mnum, &mtype, &phase,
			&rows, &ddum, CSR_I_row_indices, CSR_J_col_indices, &idum, &nRhs,
			iparm, &msglvl, &ddum, &ddum, &error, dparm);
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
			//m_error = dss_solve_real (m_handle, m_opt, &rhs[0], nRhs_l, &sol[0]);
			Solve(rhs, sol, nRhs_l);

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

	if (!keep_factors) {
		/* -------------------------------------------------------------------- */
		/* .. Termination and release of memory. */
		/* -------------------------------------------------------------------- */
		phase = -1;			/* Release internal memory. */
		int nRhs = 1;
		double ddum;
		int idum;
		pardiso (pt, &maxfct, &mnum, &mtype, &phase,
			&rows, &ddum, CSR_I_row_indices, CSR_J_col_indices, &idum, &nRhs,
			iparm, &msglvl, &ddum, &ddum, &error, dparm);
		initialized = false;
		if (MPIrank == 0) printf(".");
	}

}


void SparseSolver::SolveMat_Dense( SparseMatrix & A ) {
	SolveMat_Dense(A, A);
}

void SparseSolver::SolveMat_Dense( SparseMatrix & A_in, SparseMatrix & B_out ) {

	SEQ_VECTOR<double> rhs;
	SEQ_VECTOR<double> sol;

	eslocal job[8];
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

	eslocal m = A_in.rows;
	eslocal n = A_in.cols;
	eslocal nRhs = A_in.cols;

	rhs.resize(m * n);
	sol.resize(m * n);

	eslocal lda  = m;
	eslocal info;


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
	int nnzmax = B_out.CSR_I_row_indices[m];//-1;

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
	int iparm[64];
	int maxfct, mnum, phase, error; //, msglvl;
	/* Auxiliary variables. */
	int i;
	double ddum;			/* Double dummy */
	int idum;			/* Integer dummy. */

	/* -------------------------------------------------------------------- */
	/* .. Setup Pardiso control parameters. */
	/* -------------------------------------------------------------------- */
	for (i = 0; i < 64; i++) {
		iparm[i] = 0;
	}

	int mtype = 2;

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

	eslocal job[8];

	eslocal m		= A_in.rows;
	eslocal n		= A_in.cols;
	eslocal nRhs	= A_in.cols;
	int nRhsi       = A_in.cols;
	eslocal lda     = m;
	eslocal info;

	SEQ_VECTOR<double>  sol  (m * n, 0);

	bool clear_dense = false;

	//TODO: na konci se musi smazat !!
	if (A_in.dense_values.size() == 0) {
		A_in.ConvertCSRToDense(0);
		clear_dense = true;
	}

	SEQ_VECTOR<int> perm (A_in.dense_values.size() , 0);
	for (int ii = 0; ii < A_in.dense_values.size(); ii++)
		if (A_in.dense_values[ii] != 0.0)
			perm[ii] = 1;

	phase = 13;
	pardiso (pt, &maxfct, &mnum, &mtype, &phase,
		    &rows, CSR_V_values, CSR_I_row_indices, CSR_J_col_indices, &perm[0], &nRhsi, iparm, &msglvl, &A_in.dense_values[0], &sol[0], &error, dparm);


	if (error != 0)
	{
		printf ("\nERROR during the solution of the system : %d", error);
		exit (1);
	} else {
		initialized = true;
	}

	int x = 0;

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
	int nnzmax = B_out.CSR_I_row_indices[m];//-1; POZOR

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
	pardiso (pt, &maxfct, &mnum, &mtype, &phase,
		&rows, &ddum, CSR_I_row_indices, CSR_J_col_indices, &idum, &nRhsi,
		iparm, &msglvl, &ddum, &ddum, &error, dparm);
 
	//SEQ_VECTOR<double>().swap( A_in.dense_values );

	initialized = false;

	if (clear_dense) {
		SEQ_VECTOR<double>().swap( A_in.dense_values );
	}
}

void SparseSolver::Create_SC( SparseMatrix & B_out, int sc_size, bool isThreaded ) {

	/* Internal solver memory pointer pt, */
	/* 32-bit: int pt[64]; 64-bit: long int pt[64] */
	/* or void *pt[64] should be OK on both architectures */
	void *pt[64];

	/* Pardiso control parameters. */
	int 	iparm[64];
	double  dparm[65];
	int 	maxfct, mnum, phase, error;

	/* Auxiliary variables. */
	int 	i;
	double 	ddum;			/* Double dummy */
	int 	idum;			/* Integer dummy. */
	int 	solver;

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
	for (i = 0; i < 64; i++) {
		pt[i] = 0;
	}

	int 	mtype = 2;

//	iparm[0] = 1;		/* No solver default */
//	iparm[1] = 2;		/* Fill-in reordering from METIS */
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
//						/* Numbers of processors, value of OMP_NUM_THREADS */
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


	solver = 0; 		/* use sparse direct solver */
	pardisoinit(pt,  &mtype, &solver, iparm, dparm, &error);

    iparm[10] = 1;
    iparm[12] = 0;

 	maxfct = 1;			/* Maximum number of numerical factorizations. */
	mnum   = 1;			/* Which factorization to use. */
	//msglvl = 0;			/* Supress printing statistical information */
	error  = 0;			/* Initialize error flag */

    if (error != 0)
    {
        if (error == -10 )
           printf("No license file found \n");
        if (error == -11 )
           printf("License is expired \n");
        if (error == -12 )
           printf("Wrong username or hostname \n");
         exit(1);
    }
    else {
//        printf("[PARDISO]: License check was successful ... \n");
    }


    int nrows_S = sc_size;
    phase       = 12;
    iparm[37]   = nrows_S;

    int nb = 0; // number of righhand sides

    pardiso(pt, &maxfct, &mnum, &mtype, &phase,
               &rows,
			   CSR_V_values, CSR_I_row_indices, CSR_J_col_indices,
			   &idum, &nb,
               iparm, &msglvl, &ddum, &ddum, &error,  dparm);

    if (error != 0)
    {
        printf("\nERROR during symbolic factorization: %d", error);
        exit(1);
    } else {
    	initialized = true;
    }
//    printf("\nReordering completed ...\n");
//    printf("Number of nonzeros in factors  = %d\n", iparm[17]);
//    printf("Number of factorization MFLOPS = %d\n", iparm[18]);
//    printf("Number of nonzeros is   S      = %d\n", iparm[38]);

    /* -------------------------------------------------------------------- */
    /* ..  allocate memory for the Schur-complement and copy it there.      */
    /* -------------------------------------------------------------------- */
    int nonzeros_S = iparm[38];

    B_out.CSR_I_row_indices.resize(nrows_S+1);
    B_out.CSR_J_col_indices.resize(nonzeros_S);
    B_out.CSR_V_values.resize(nonzeros_S);
    B_out.cols = nrows_S;
    B_out.rows = nrows_S;
    B_out.nnz  = nonzeros_S;
    B_out.type = 'S';

    pardiso_get_schur(pt, &maxfct, &mnum, &mtype, &B_out.CSR_V_values[0], &B_out.CSR_I_row_indices[0], &B_out.CSR_J_col_indices[0]);

    phase = -1;                 /* Release internal memory. */

    pardiso(pt, &maxfct, &mnum, &mtype, &phase,
             &rows, &ddum, CSR_I_row_indices, CSR_J_col_indices, &idum, &idum,
             iparm, &msglvl, &ddum, &ddum, &error,  dparm);

    initialized = false;

}

void SparseSolver::Create_SC_w_Mat( SparseMatrix & K_in, SparseMatrix & B_in, SparseMatrix & SC_out,
								    bool isThreaded, int generate_symmetric_sc_1_generate_general_sc_0 ) {

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
	int 	iparm[64];
	double  dparm[65];
	int 	maxfct, mnum, phase, error;

	/* Auxiliary variables. */
	int 	i;
	double 	ddum;			/* Double dummy */
	int 	idum;			/* Integer dummy. */
	int 	solver;

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
	for (i = 0; i < 64; i++) {
		pt[i] = 0;
	}

	int 	mtype = 2;

//	iparm[0] = 1;		/* No solver default */
//	iparm[1] = 2;		/* Fill-in reordering from METIS */
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
//						/* Numbers of processors, value of OMP_NUM_THREADS */
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


	solver = 0; 		/* use sparse direct solver */
	pardisoinit(pt,  &mtype, &solver, iparm, dparm, &error);

    iparm[10] = 1;
    iparm[12] = 0;

 	maxfct = 1;			/* Maximum number of numerical factorizations. */
	mnum   = 1;			/* Which factorization to use. */
	//msglvl = 0;		/* Supress printing statistical information */
	error  = 0;			/* Initialize error flag */

    if (error != 0)
    {
        if (error == -10 )
           printf("No license file found \n");
        if (error == -11 )
           printf("License is expired \n");
        if (error == -12 )
           printf("Wrong username or hostname \n");
         exit(1);
    }
    else {
//        printf("[PARDISO]: License check was successful ... \n");
    }


    int nrows_S = B_in.cols;
    phase       = 12;
    iparm[37]   = nrows_S;

    int nb = 0; // number of righhand sides

    pardiso(pt, &maxfct, &mnum, &mtype, &phase,
               &K_sc1.rows,
			   &K_sc1.CSR_V_values[0], &K_sc1.CSR_I_row_indices[0], &K_sc1.CSR_J_col_indices[0],
			   &idum, &nb,
               iparm, &msglvl, &ddum, &ddum, &error,  dparm);

    if (error != 0)
    {
        printf("\nERROR during symbolic factorization: %d", error);
        exit(1);
    } else {
    	initialized = true;
    }

    /* -------------------------------------------------------------------- */
    /* ..  allocate memory for the Schur-complement and copy it there.      */
    /* -------------------------------------------------------------------- */
    int nonzeros_S = iparm[38];

    SC_out.CSR_I_row_indices.resize(nrows_S+1);
    SC_out.CSR_J_col_indices.resize(nonzeros_S);
    SC_out.CSR_V_values.resize(nonzeros_S);
    SC_out.cols = nrows_S;
    SC_out.rows = nrows_S;
    SC_out.nnz  = nonzeros_S;
    SC_out.type = 'S';

    pardiso_get_schur(pt, &maxfct, &mnum, &mtype, &SC_out.CSR_V_values[0], &SC_out.CSR_I_row_indices[0], &SC_out.CSR_J_col_indices[0]);

    phase = -1;                 /* Release internal memory. */

    pardiso(pt, &maxfct, &mnum, &mtype, &phase,
    		&K_sc1.rows, &ddum, &K_sc1.CSR_I_row_indices[0], &K_sc1.CSR_J_col_indices[0], &idum, &idum,
             iparm, &msglvl, &ddum, &ddum, &error,  dparm);

    initialized = false;

    // Finalize shape of the SC
    if (generate_symmetric_sc_1_generate_general_sc_0 == 0) {

		SC_out.type = 'G';

		SparseMatrix SC_tmp;
		SC_tmp = SC_out;
		SC_tmp.SetDiagonalOfSymmetricMatrix(0.0);
		SC_tmp.MatTranspose();

		SC_out.MatAddInPlace(SC_tmp,'N',1.0);
		SC_out.MatScale(-1.0);
		SC_out.ConvertCSRToDense(1);

    } else {
		SC_out.MatScale(-1.0);
		SC_out.ConvertCSRToDense(1);
		SC_out.RemoveLowerDense();
    }


}

void SparseSolver::Create_non_sym_SC_w_Mat( SparseMatrix & K_in, SparseMatrix & B1_in, SparseMatrix & B0_in, SparseMatrix & SC_out, bool isThreaded, int generate_symmetric_sc_1_generate_general_sc_0 ) {

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
	int 	iparm[64];
	double  dparm[65];
	int 	maxfct, mnum, phase, error;

	/* Auxiliary variables. */
	int 	i;
	double 	ddum;			/* Double dummy */
	int 	idum;			/* Integer dummy. */
	int 	solver;

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
	for (i = 0; i < 64; i++) {
		pt[i] = 0;
	}

	int 	mtype = 11; // real unsymetric

//	iparm[0] = 1;		/* No solver default */
//	iparm[1] = 2;		/* Fill-in reordering from METIS */
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
//						/* Numbers of processors, value of OMP_NUM_THREADS */
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


	solver = 0; 		/* use sparse direct solver */
	pardisoinit(pt,  &mtype, &solver, iparm, dparm, &error);

    iparm[10] = 1;
    iparm[12] = 0;

 	maxfct = 1;			/* Maximum number of numerical factorizations. */
	mnum   = 1;			/* Which factorization to use. */
	//msglvl = 1;		/* Supress printing statistical information */
	error  = 0;			/* Initialize error flag */

    if (error != 0)
    {
        if (error == -10 )
           printf("No license file found \n");
        if (error == -11 )
           printf("License is expired \n");
        if (error == -12 )
           printf("Wrong username or hostname \n");
         exit(1);
    }
    else {
//        printf("[PARDISO]: License check was successful ... \n");
    }


    int nrows_S = B1_in.cols;
    phase       = 12;
    iparm[37]   = nrows_S;

    int nb = 0; // number of righhand sides

    pardiso(pt, &maxfct, &mnum, &mtype, &phase,
               &K_sc1.rows,
			   &K_sc1.CSR_V_values[0], &K_sc1.CSR_I_row_indices[0], &K_sc1.CSR_J_col_indices[0],
			   &idum, &nb,
               iparm, &msglvl, &ddum, &ddum, &error,  dparm);

    if (error != 0)
    {
        printf("\nERROR during symbolic factorization: %d", error);
        exit(1);
    } else {
    	initialized = true;
    }

    /* -------------------------------------------------------------------- */
    /* ..  allocate memory for the Schur-complement and copy it there.      */
    /* -------------------------------------------------------------------- */
    int nonzeros_S = iparm[38];

    SC_out.CSR_I_row_indices.resize(nrows_S+1);
    SC_out.CSR_J_col_indices.resize(nonzeros_S);
    SC_out.CSR_V_values.resize(nonzeros_S);
    SC_out.cols = nrows_S;
    SC_out.rows = nrows_S;
    SC_out.nnz  = nonzeros_S;
    SC_out.type = 'G';

    pardiso_get_schur(pt, &maxfct, &mnum, &mtype, &SC_out.CSR_V_values[0], &SC_out.CSR_I_row_indices[0], &SC_out.CSR_J_col_indices[0]);

    phase = -1;                 /* Release internal memory. */

    pardiso(pt, &maxfct, &mnum, &mtype, &phase,
    		&K_sc1.rows, &ddum, &K_sc1.CSR_I_row_indices[0], &K_sc1.CSR_J_col_indices[0], &idum, &idum,
             iparm, &msglvl, &ddum, &ddum, &error,  dparm);

    initialized = false;

    // Finilize shape of the SC

    //SpyText(SC_out);

    SC_out.ConvertToCOO(1);
    SC_out.rows = B0_in.cols;

    vector <int> 	I;
    vector <int> 	J;
    vector <double> V;

    for (int i = 0; i < SC_out.I_row_indices.size(); i++) {
    	if ( SC_out.I_row_indices[i] <= B0_in.cols) {
    		I.push_back(SC_out.I_row_indices[i]);
    		J.push_back(SC_out.J_col_indices[i]);
    		V.push_back(SC_out.V_values[i]);
    	}
    }

    SC_out.I_row_indices = I;
    SC_out.J_col_indices = J;
    SC_out.V_values      = V;

    SC_out.nnz = V.size();

    SC_out.ConvertToCSRwithSort(1);

    //SpyText(SC_out);

    SC_out.MatScale(-1.0);

}

