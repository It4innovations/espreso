
#ifndef SOLVER_SPECIFIC_CPU_SOLVERDISSECTION_H_
#define SOLVER_SPECIFIC_CPU_SOLVERDISSECTION_H_

#include "densesolvers.h"

#include "sparsesolver.h"

#include "mkl_pardiso.h"

// #include "mkl_types.h"
// #undef MKL_Complex16

template<typename T, typename U, typename W, typename Z, typename X, typename Y>
class DissectionSolver;


namespace espreso {

class SparseSolverDissection: public SparseSolver
{

public:
	SparseSolverDissection();
	~SparseSolverDissection();


	//Members
	void ImportMatrix(SparseMatrix & A);
	void ImportMatrix_fl(SparseMatrix & A);

	void ImportMatrix_wo_Copy(SparseMatrix & A);
	void ImportMatrix_wo_Copy_fl(SparseMatrix & A);

	void ExportMatrix(espreso::SparseMatrix & A);

	int Factorization(const std::string &str);
	void Clear();
	void SetThreaded();


	void Solve( std::vector <double> & rhs, std::vector <double> & sol, MKL_INT rhs_start_index, MKL_INT sol_start_index);
	void Solve( std::vector <double> & rhs, std::vector <double> & sol, MKL_INT n_rhs );
	void Solve( std::vector <double> & rhs_sol);

	void SolveMat_Dense( SparseMatrix & A_in_out );
	void SolveMat_Dense( SparseMatrix & A_in, SparseMatrix & B_out );

	void SolveMatF( SparseMatrix & A_in, SparseMatrix & B_out, bool isThreaded );

	void SolveMat_Sparse( SparseMatrix & A );
	void SolveMat_Sparse( SparseMatrix & A_in, SparseMatrix & B_out );
	void SolveMat_Sparse( SparseMatrix & A_in, SparseMatrix & B_out, char T_for_input_matrix_is_transposed_N_input_matrix_is_NOT_transposed );

	void Create_SC( SparseMatrix & B_out, MKL_INT sc_size, bool isThreaded );
	void Create_SC_w_Mat( SparseMatrix & K_in, SparseMatrix & B_in, SparseMatrix & SC_out, bool isThreaded, MKL_INT generate_symmetric_sc_1_generate_general_sc_0 );
	void Create_non_sym_SC_w_Mat( SparseMatrix & K_in, SparseMatrix & B1_in, SparseMatrix & B0_in, SparseMatrix & SC_out, bool isThreaded, MKL_INT generate_symmetric_sc_1_generate_general_sc_0 );

	void SolveCG(SparseMatrix & A_in, SEQ_VECTOR <double> & rhs_in, SEQ_VECTOR <double> & sol, SEQ_VECTOR <double> & initial_guess);
	void SolveCG(SparseMatrix & A_in, std::vector <double> & rhs, std::vector <double> & sol);
	void SolveCG(SparseMatrix & A_in, std::vector <double> & rhs_sol);

//	void GetKernelVectors(SEQ_VECTOR <double> & kern_vec, esint & kern_dim);
//	void GetTransKernelVectors(SEQ_VECTOR <double> & kern_t_vec, esint & kern_dim);
	void GetKernel(SparseMatrix &R, SparseMatrix &R2);

	void SaveMatrixInCSR(string filename);

	DissectionSolver<double, double, double, double, double, double> * dslv;
	FILE *fp;
	int called;
	double eps_pivot;
	int scaling;
	bool kernel_detection_all;
	bool diss_verbose;
	bool projection;

	bool use_dense_solver;
	DenseSolverCPU dense_solver;

	bool initialized;
	bool keep_factors;
	bool import_with_copy;
	bool USE_FLOAT;
	int  MPIrank;
	int num_procs;

	MKL_INT rows;
	MKL_INT cols;
	MKL_INT nnz;

	bool is_whole;
	bool is_sym;
	bool upper_flag;
	int decomposer;
	int nb_levels;
	int min_nodes;


	MKL_INT		* I_row_indices;
	MKL_INT		* J_col_indices;
	double		* V_values;
	float		* V_values_fl;

	MKL_INT		I_row_indices_size;
	MKL_INT		J_col_indices_size;
	MKL_INT		V_values_size;


	MKL_INT		* CSR_I_row_indices;
	MKL_INT		* CSR_J_col_indices;
	double		* CSR_V_values;
	float 		* CSR_V_values_fl;

	MKL_INT		CSR_I_row_indices_size;
	MKL_INT		CSR_J_col_indices_size;
	MKL_INT		CSR_V_values_size;
	MKL_INT 	CSR_V_values_fl_size;

	// *** Pardiso Solver Variables

	MKL_INT mtype;		/* Real symmetric matrix */

	/* Internal solver memory pointer pt, */
	/* 32-bit: int pt[64]; 64-bit: long int pt[64] */
	/* or void *pt[64] should be OK on both architectures */
	// void *pt[64];
	//int * pt;


	/* Pardiso control parameters. */
	MKL_INT iparm[65];
	double  dparm[65];

	// MKL_INT maxfct, mnum, phase, error, msglvl;
	MKL_INT error, msglvl;
	// ***

	MKL_INT m_nRhs;
	MKL_INT m_factorized;
	MKL_INT m_Kplus_size;
	// END - MKL DSS Solver Variables

	// Matrices

	// for in-place solve
	std::vector <double> tmp_sol;
	std::vector <float> tmp_sol_fl1;
	std::vector <float> tmp_sol_fl2;

};

}


#endif /* SOLVER_SPECIFIC_CPU_SOLVERDISSECTION_H_ */
