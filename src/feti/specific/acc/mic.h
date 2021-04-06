
#ifndef SOLVER_SPARSE_MIC_MKL_H_
#define SOLVER_SPARSE_MIC_MKL_H_

#include "sparsesolver.h"
#include "mkl_pardiso.h"
#include "DenseMatrixPack.h"
#include "SparseMatrixPack.h"

namespace espreso {

class SparseSolverMIC
{

public:
	SparseSolverMIC();
	~SparseSolverMIC();


	//Members
	void ImportMatrices(SparseMatrix ** A, MKL_INT nMatrices, esint mic = 0);
	void ImportMatrices_fl(SparseMatrix ** A, esint nMatrices, esint mic = 0);

	void ImportMatrices_wo_Copy(SparseMatrix ** A, esint nMatrices, esint mic = 0 );

	void Factorization(const std::string &str);
    void Factorization(const std::string &str, SparseMatrixPack & factors_out);
	void Clear();
	void SetThreaded() {};


	void Solve( std::vector <double> & rhs, std::vector <double> & sol, MKL_INT rhs_start_index, MKL_INT sol_start_index) {};
	void Solve( std::vector <double> & rhs, std::vector <double> & sol, MKL_INT n_rhs ) {};
	void Solve( std::vector <double> ** rhs_sol);

	void SolveMat_Dense( SparseMatrix & A_in_out ) {};
	void SolveMat_Dense( SparseMatrix & A_in, SparseMatrix & B_out ) {};

	void SolveMatF( SparseMatrix & A_in, SparseMatrix & B_out, bool isThreaded ) {};

	void SolveMat_Sparse( SparseMatrix & A ) {};
	void SolveMat_Sparse( SparseMatrix & A_in, SparseMatrix & B_out ) {};
	void SolveMat_Sparse( SparseMatrix & A_in, SparseMatrix & B_out, char T_for_input_matrix_is_transposed_N_input_matrix_is_NOT_transposed ) {};

	void Create_SC( SparseMatrix & B_out, MKL_INT sc_size, bool isThreaded ) {};
    void Create_SC( DenseMatrixPack & SC_out, MKL_INT* sc_sizes, esint generate_symmetric_sc_1_generate_general_sc_0 );
	void Create_SC_w_Mat( SparseMatrix & K_in, SparseMatrix & B_in, SparseMatrix & SC_out, bool isThreaded, MKL_INT generate_symmetric_sc_1_generate_general_sc_0 ) {};
    void Create_SC_w_Mat( SparseMatrix ** K_in, SparseMatrix ** B_in, DenseMatrixPack & SC_out, esint nMatrices, esint generate_symmetric_sc_1_generate_general_sc_0, esint device  );
	void Create_non_sym_SC_w_Mat( SparseMatrix & K_in, SparseMatrix & B1_in, SparseMatrix & B0_in, SparseMatrix & SC_out, bool isThreaded, MKL_INT generate_symmetric_sc_1_generate_general_sc_0 ) {};

	void SolveCG(SparseMatrix & A_in, std::vector <double> & rhs, std::vector <double> & sol) {};
	void SolveCG(SparseMatrix & A_in, std::vector <double> & rhs_sol) {};

    esint device; // which mic to use for offload
    bool isOffloaded;

    esint nMatrices;
	bool initialized;
	bool keep_factors;
	bool import_with_copy;
	bool USE_FLOAT;
	esint  MPIrank;

	MKL_INT* rows;
	MKL_INT* cols;
	MKL_INT* nnz;



	MKL_INT		** I_row_indices;
	MKL_INT		** J_col_indices;
	double		** V_values;
	float		** V_values_fl;

	MKL_INT		* I_row_indices_size;
	MKL_INT		* J_col_indices_size;
	MKL_INT		* V_values_size;


	MKL_INT		** CSR_I_row_indices;
	MKL_INT		** CSR_J_col_indices;
	double		** CSR_V_values;
	float 		** CSR_V_values_fl;

	MKL_INT		* CSR_I_row_indices_size;
	MKL_INT		* CSR_J_col_indices_size;
	MKL_INT		* CSR_V_values_size;
	MKL_INT 	* CSR_V_values_fl_size;

	// *** Pardiso Solver Variables

	MKL_INT mtype;		/* Real symmetric matrix */

	/* Internal solver memory poesinter pt, */
	/* 32-bit: esint pt[64]; 64-bit: long esint pt[64] */
	/* or void *pt[64] should be OK on both architectures */
	void ***pt;
	//esint * pt;


	/* Pardiso control parameters. */
	MKL_INT **iparm;
    double  **dparm;
	MKL_INT ** perm;

    // common parameters
	MKL_INT maxfct, mnum, phase, msglvl;
	// ***
    MKL_INT *error;


	MKL_INT m_nRhs;
	MKL_INT m_factorized;
	MKL_INT *m_Kplus_size;
	// END - MKL DSS Solver Variables

	// Matrices
	//SparseMatrix m_A;

    SparseMatrix ** A;

	// for in-place solve
	// std::vector <double> tmp_sol;
	float ** tmp_sol_fl1;
	float ** tmp_sol_fl2;
    double ** tmp_sol_d1;
    double ** tmp_sol_d2;

    // array of arrays of all vectors
    double * vectors;
    double * vectors_out;
};

}

#endif
