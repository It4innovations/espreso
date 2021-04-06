/** @file SparseSolverMumps.cpp
	@brief Solver using MUMPS library.

	@author Martin Beseda
	@date 2015-

	@todo Dealloc MUMPS
*/

/************************************************************************************************
*   INCLUDE LIBRARIES   *
************************/
#include "mkl_pardiso.h"
#include "solvermumps.h"

#include "basis/utilities/utils.h"

using namespace espreso;

/************************************************************************************************
*   DEFINE CONSTANTS   *
***********************/
#define JOB_INIT -1			   /**< id%JOB value for MUMPS instance - initialize instance. */
#define JOB_END -2			   /**< id%JOB value for MUMPS instance - terminate instance. */
#define USE_COMM_WORLD -987654 /**< info::mpi::MPICommunicator for MUMPS instance.  */

/************************************************************************************************
*   DEFINE MACROS   *
********************/
#define ICNTL(I) icntl[(I)-1]

/************************************************************************************************
*   SparseSolver METHODS   *
***************************/
/** @brief Controls, if the variable SOLVER_NUM_THREADS is set.
 *
 */
void SparseSolverMUMPS::SetThreaded() {

	/* Numbers of processors, value of OMP_NUM_THREADS */
	int num_procs = Esutils::getEnv<int>("SOLVER_NUM_THREADS");

}

/** @brief Computes Schur complement.
 *
 * @param[out] B_out      SparseMatrix, where Schur complement is saved.
 * @param[in]  sc_size    Schur complement size.
 */
void SparseSolverMUMPS::Create_SC( SparseMatrix & B_out, int sc_size, bool isThreaded ) {
	/* Note, that Schur complement is stored BY COLUMNS,
	   storing by rows is obsolete in MUMPS. */

	id.job			 = 1; // analysis

	id.ICNTL(19)	 = 2; // 2 - returns only lower triangular matrix if symmetric, 3 - returns whole matrix if symmetric
	id.size_schur	 = sc_size; // size (order) of Schur complement
	id.listvar_schur = new int[sc_size]; // values of matrix 'id.a' which are used for computing Schur c.

	// Init 'id.listvar_schur' with indices of columns from 'id.a' used for computing Schur complement
    int schur_base_index = id.n - id.size_schur;
	for(int i = 1; i <= sc_size; i++) {
		id.listvar_schur[i] = schur_base_index + i;
	};

	id.nprow		 = 1;
	id.npcol		 = 1;
	id.mblock		 = 10;//100;
	id.nblock		 = 10;//100;

	dmumps_c(&id);

	id.job			 = 2; // factorization

	id.schur_lld	 = sc_size;

	// INIT output matrix 'B_out'
	B_out.rows = sc_size;
	B_out.cols = sc_size;
	B_out.type = 'S'; // 'S' - symmetric, 'G' - general - WARNING: must correspond with id.ICNTL(19) value!
	B_out.dense_values.resize(sc_size * sc_size);

	id.schur		 = &B_out.dense_values[0]; // assign array for output

	dmumps_c(&id);
}

void SparseSolverMUMPS::Create_SC_w_Mat( SparseMatrix & K_in, SparseMatrix & B_in, SparseMatrix & SC_out, bool isThreaded, MKL_INT generate_symmetric_sc_1_generate_general_sc_0 )
{
	ESINFO(ERROR) << "Not implemented in MUMPS";
}
/** @brief Constructor.
 *
 */
SparseSolverMUMPS::SparseSolverMUMPS(){

	I_row_indices_size = 0;
	J_col_indices_size = 0;
	V_values_size = 0;

	id.ICNTL(2) = 0;
	id.ICNTL(3) = 0;

	id.job          = JOB_INIT;
	id.par          = 1;
	id.sym          = 1; // 0 - unsymmetric matrix, 1 - positive definite matrix, 2 - generally symetric matrix
	id.comm_fortran = (MUMPS_INT) MPI_Comm_c2f(MPI_COMM_SELF);

	dmumps_c(&id);

	maxfct = 1;				/* Maximum number of numerical factorizations. */
	//mnum = 1;				/* Which factorization to use. */
#ifdef DEBUG
	msglvl = 1;				/* Print statistical information in file */
#else
	msglvl = 0;
#endif
	error = 0;				/* Initialize error flag */

	m_nRhs		 = 1;
	m_factorized = 0;
}

/** @brief Desctructor.
 *
 * Deallocates I_row_indices, J_col_indices, V_values
 */
SparseSolverMUMPS::~~SparseSolverMUMPS() {

	// MUMPS instance termination TODO dodelat - MUMPS se musi dealokovat
	//id.job=JOB_END;
	//dmumps_c(&id);

	if (I_row_indices_size > 0)		delete [] I_row_indices;
	if (J_col_indices_size > 0)		delete [] J_col_indices;
	if (V_values_size > 0)			delete [] V_values;

	I_row_indices_size = 0;
	J_col_indices_size = 0;
	V_values_size      = 0;
}

/** @brief Deallocs matrix, like the desctructor.
 *
 * @see SparseSolverMUMPS::~SparseSolver()
 */
void SparseSolverMUMPS::Clear() {
	/* -------------------------------------------------------------------- */
	/* .. Termination and release of memory. */
	/* -------------------------------------------------------------------- */

	//TODO dodelat MUMPS

	if (I_row_indices_size > 0)     delete [] I_row_indices;
	if (J_col_indices_size > 0)		delete [] J_col_indices;
	if (V_values_size > 0)			delete [] V_values;

	I_row_indices_size = 0;
	J_col_indices_size = 0;
	V_values_size      = 0;
}

/** @brief Imports matrix to the solver.
 *
 * @param[in] A SparseMatrix to be imported.
 */
void SparseSolverMUMPS::ImportMatrix(SparseMatrix & A_in) {

	SparseMatrix A;
	A = A_in;

	rows	= A.rows;
	cols	= A.cols;
	nnz		= A.nnz;

	m_Kplus_size = A.rows;

	A.ConvertToCOO(0); // 0 - keep CSR, 1 - remove CSR // TODO poresit, jestli mazat

	I_row_indices_size = A.I_row_indices.size();
	J_col_indices_size = A.J_col_indices.size();
	V_values_size	   = A.V_values.size();

	I_row_indices = new int[I_row_indices_size];
	J_col_indices = new int[J_col_indices_size];
	V_values	  = new double[V_values_size];

	copy(A.I_row_indices.begin(), A.I_row_indices.end(), I_row_indices);
	copy(A.J_col_indices.begin(), A.J_col_indices.end(), J_col_indices);
	copy(A.V_values     .begin(), A.V_values     .end(), V_values);

	// Init MUMPS instance
	id.n   = A.rows; // matrix order
	id.nz  = A.nnz; // number of nonzero elements
	id.irn = I_row_indices; // horizontal coordinates of nonzero elements
	id.jcn = J_col_indices; // vertical coordinates of nonzero elements
	id.a   = V_values; // the matrix itself

}

/** @brief Does factorization of imported matrix.
 *
 */
void SparseSolverMUMPS::Factorization() {
	//phase = 11;

	//SparseMatrix B;
	//Create_SC(B, 5, false);

	id.job = 1;
	dmumps_c(&id); // Analysis

	id.job = 2;
	dmumps_c(&id); // Factorization

#ifdef DEBUG
	ESINFO(PROGRESS3) << "Factorization completed ... ";
#endif

	m_factorized = 1;
}

/** @brief Solves system of equations with one right-hand side.
 *
 * @param[in,out] rhs_sol Vector containing right-hand side vector as input
 *                        and the solutioni of the system as output.
 */
void SparseSolverMUMPS::Solve(SEQ_VECTOR <double> & rhs_sol) {
	id.job = 3; // set MUMPS to solve the system

	id.rhs = &rhs_sol[0]; // init MUMPS with right-hand side vector

	dmumps_c(&id); // solve the system

#ifdef DEBUG
	printf ("\nSolve completed ... ");
	printf ("\nThe solution of the system is: ");
	for (int i = 0; i < n; i++)
	{
		printf ("\n x [%d] = % f", i, x[i]);
	}
	printf ("\n");
#endif
}

/** @brief Solves system of equations with multiple right-hand sides.
 *
 * @param[out]    sol   Vector with solution of the system.
 * @param[in]     rhs   Multiple right-hand sides of the system.
 * @param[in,out] n_rhs Number of right-hand side vectors.
 */
void SparseSolverMUMPS::Solve(SEQ_VECTOR <double> & rhs, SEQ_VECTOR <double> & sol, int n_rhs) {
	id.job		   = 3; // set MUMPS to solve the system

	int rhs_size   = rhs.size() / n_rhs; // size of one InitialCondition

    sol = rhs;

	id.ICNTL(20) = 0; // dense InitialCondition
	id.nrhs		 = n_rhs; // number of InitialCondition vectors
	id.lrhs		 = rhs_size; // size of one InitialCondition vector
	id.rhs		 = &sol[0]; // init MUMPS instance with vector of InitialCondition ('rhs' param)

	dmumps_c(&id); // solve the system

#ifdef DEBUG
	printf ("\nSolve completed ... ");
	printf ("\nThe solution of the system is: ");
	for (int i = 0; i < n; i++)
	{
		printf ("\n x [%d] = % f", i, x[i]);
	}
	printf ("\n");
#endif
}

/** @brief Solves the system of equations, with multiple right-hand sides,
 *         where the user can choose the starting index of InitialCondition.
 *
 * TODO poresit nejasnosti s indexy a dopsat jejich popis
 * @param[out] sol Solution of the system
 * @param[in]  rhs Right-hand sides of the system.
 * @param[in]  rhs_start_index
 * @param[in]  sol_start_index
 */
void SparseSolverMUMPS::Solve(SEQ_VECTOR <double> & rhs, SEQ_VECTOR <double> & sol, int rhs_start_index, int sol_start_index) {
    sol = rhs;

	id.job = 3; // set MUMPS to solve the system
	id.rhs = &sol[rhs_start_index];

	dmumps_c(&id);

#ifdef DEBUG
	printf ("\nSolve completed ... ");
	printf ("\nThe solution of the system is: ");
	for (int i = 0; i < n; i++)
	{
		printf ("\n x [%d] = % f", i, x[i]);
	}
	printf ("\n");
#endif
}

/** @brief Solves the system of equations where solutions are returned as a SparseMatrix.
 *
 * @param[in,out] A SparseMatrix which contains multiple right-hand sides (column-wise)
 *                  for input and solutions (column-wise) for output.
 */
void SparseSolverMUMPS::SolveMat_Sparse( SparseMatrix & A) {
	SolveMat_Sparse(A, A);
};

/** @brief Solves the system of equations where solutions are returned as a SparseMatrix.
 *
 * @param[out] B_out SparseMatrix containing multiple solutions (column-wise).
 * @param[in]  A_in  SparseMatrix containing multiple right-hand sides of the system, for
 *                   input is transposed!
 */
void SparseSolverMUMPS::SolveMat_Sparse( SparseMatrix & A_in, SparseMatrix & B_out) {
	SolveMat_Sparse(A_in, B_out, 'T');
};

/** @brief Solves the system of equations where solutions are returned as a SparseMatrix.
 *
 * @param[out] B_out SparseMatrix containing all solutions (column-wise).
 * @param[in]  A_in  SparseMatrix containing all right-hand sides (column-wise).
 * @param[in]  T_for_input_matrix_is_transposed_N_input_matrix_is_NOT_transposed The switch which allows to transpose the input matrix before solving the system.
 */
void SparseSolverMUMPS::SolveMat_Sparse( SparseMatrix & A_in, SparseMatrix & B_out, char T_for_input_matrix_is_transposed_N_input_matrix_is_NOT_transposed ) {
	// TODO popsat system prevodu CSR na CSC pomoci maticove transpozice
	SparseMatrix tmpM;

    if (T_for_input_matrix_is_transposed_N_input_matrix_is_NOT_transposed == 'T') {
        //A_in.MatTranspose(tmpM);
		tmpM = A_in;
	} else {
        //tmpM = A_in;
		A_in.MatTranspose(tmpM); // because of "transfer" from CSR to CSC
	}

    /* Matrix 'A_in' processed in CSC (Compressed Sparse Column) format */
	id.ICNTL(20)   = 1; // 2 - exploit sparsity to improve solution phase, 3 - not exploit sparsity, 1 - decide automatically
	id.irhs_ptr	   = &tmpM.CSR_I_row_indices[0]; // "pointers" to the columns
	id.irhs_sparse = &tmpM.CSR_J_col_indices[0]; // row indices
	id.rhs_sparse  = &tmpM.CSR_V_values[0]; // non-zero values
	id.nz_rhs	   = tmpM.nnz; // total amount of non-zero values in InitialCondition
	id.nrhs		   = tmpM.cols; // number of InitialCondition vectors
    id.lrhs        = tmpM.rows; // TODO zjistit, jestli je nutne

	unsigned int id_rhs_len = A_in.cols * A_in.rows;
	id.rhs		   = new double[id_rhs_len]; // output in dense format

	id.job		   = 3; // solve the system
	dmumps_c(&id);

	// Create and init output matrix 'B_out'
	B_out.rows = A_in.rows;
	B_out.cols = A_in.cols;
	B_out.type = 'G';
	B_out.dense_values.assign(id.rhs, id.rhs+id_rhs_len);

	B_out.ConvertDenseToCSR(1); // 0 - keep dense values, 1 - clear dense values

	tmpM.Clear(); // dealloc 'tmpM' matrix
}

/** @brief Solves the system of equations where solutions are returned as a dense matrix.
 *
 * @param[in,out] A SparseMatrix containing multiple right-hand sides for input and
 *                  multiple solutions in 'dense_values' attribute for output.
 */
void SparseSolverMUMPS::SolveMat_Dense( SparseMatrix & A ) {
	SolveMat_Dense(A, A);
}

/** @brief Solves the system of equations where solutions are returned as a dense matrix.
 *
 * @param[out] B_out SparseMatrix containing multiple solutions in 'dense_values' attribute.
 * @param[in]  A_in  SparseMatrix containing multiple right-hand sides column-wise.
 */
void SparseSolverMUMPS::SolveMat_Dense( SparseMatrix & A_in, SparseMatrix & B_out ) {
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


	// Convert input matrix (InitialCondition) to dense format

	//void mkl_ddnscsr (
	//	MKL_INT *job,
	//	MKL_INT *m, MKL_INT *n,
	//	double *Adns, MKL_INT *lda,
	//	double *Acsr, MKL_INT *AJ, MKL_INT *AI,
	//	MKL_INT *info);

//TODO sparse InitialCondition misto dense - i v SolveMatSparse
	mkl_ddnscsr (
		job,
		&m, &n,
		&rhs[0], &lda,
		&A_in.CSR_V_values[0], &A_in.CSR_J_col_indices[0], &A_in.CSR_I_row_indices[0],
		&info);

	// Solve with multiple right hand sides
	//	m_error = dss_solve_real (m_handle, m_opt, &rhs[0], nRhs, &sol[0]);

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

void SparseSolverMUMPS::SolveMatF( SparseMatrix & A_in, SparseMatrix & B_out, bool isThreaded) {

 ESINFO(ERROR) << "Not Implemented in MUMPS";

}




