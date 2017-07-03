
#include "../generic/utils.h"
#include "../../basis/matrices/matrixtype.h"
#include "../../configuration/environment.h"
#include "../../configuration/solver/espreso.h"

//class SparseSolverCPU;

#ifdef CUDA
	#include <cuda_runtime.h>
	#include <cublas_v2.h>
#endif

using std::vector;
using std::string;
using std::ifstream;
using std::stringstream;
using std::endl;
using std::getline;
using std::fill;


#pragma once

namespace espreso {

template <typename TIndices> class SparseCSRMatrix;
template <typename TIndices> class SparseIJVMatrix;

class SparseMatrix
{

public:

	//Constructors
	SparseMatrix(char matrix_type_G_for_general_S_for_symmetric, string filename);
	SparseMatrix(char matrix_type_G_for_general_S_for_symmetric);
	SparseMatrix();
	SparseMatrix( const SparseMatrix &A_in);
	SparseMatrix( const SparseCSRMatrix<eslocal> &A_in, char type_in );
	SparseMatrix( const SparseIJVMatrix<eslocal> &A_in, char type_in );

	friend std::ostream& operator<<(std::ostream& os, const SparseMatrix &m);


	//Destructor
	~SparseMatrix();

	//assignment operator
	SparseMatrix& operator= ( const SparseMatrix &A_in);
	SparseMatrix& operator= ( const SparseCSRMatrix<eslocal> &A_in );
	SparseMatrix& operator= ( const SparseIJVMatrix<eslocal> &A_in );

	void swap ( SparseMatrix &A_in);

	// Variables
	eslocal  rows;		// number of rows
	eslocal  cols;		// number of columns
	eslocal  nnz;		// number of non zero elements
	char type;		// 'G' for general or 'S' for symmetric
	MatrixType mtype;
	char uplo; 		// 'L' for lower or 'U' for upper
	eslocal extern_lda; // leading dimension of the array holding 2 SCs in the device memory

	// Sparse COO data
	SEQ_VECTOR <eslocal>	I_row_indices;
	SEQ_VECTOR <eslocal>	J_col_indices;
	SEQ_VECTOR <double> 	V_values;

	// Sparse CSR data
	SEQ_VECTOR <eslocal>	CSR_I_row_indices;
	SEQ_VECTOR <eslocal>	CSR_J_col_indices;
	SEQ_VECTOR <double> 	CSR_V_values;

	// Dense data
	SEQ_VECTOR <double> 	dense_values;
	SEQ_VECTOR <float>  	dense_values_fl;
	SEQ_VECTOR <eslocal>    ipiv;

	SEQ_VECTOR <float> 		vec_fl_in;
	SEQ_VECTOR <float> 		vec_fl_out;
	bool					USE_FLOAT;

	// CUDA
	double * d_dense_values;
	double * d_x_in;
	double * d_y_out;

	float  * d_dense_values_fl;
	float  * d_x_in_fl;
	float  * d_y_out_fl;

#ifdef CUDA
	cublasHandle_t handle;
	cudaStream_t stream;

	void SetCUDA_Stream(cudaStream_t & in_stream);
	void ClearCUDA_Stream();
	const char * _cudaGetErrorEnum(cublasStatus_t error);
#endif


	//Members
	eslocal  SaveMatrixInCOO(string filename);
	eslocal  SaveMatrixBinInCOO(string filename);
	eslocal  SaveMatrixBinInCSR(string filename);

	//void LoadMatrixInCOO(string filename, char matrix_type_G_for_general_S_for_symmetric);
	//void LoadMatrix(string filename, char matrix_type_G_for_general_S_for_symmetric);
	eslocal  LoadMatrixBin(string filename, char matrix_type_G_for_general_S_for_symmetric);
	eslocal  LoadMatrixBin(string filename, char matrix_type_G_for_general_S_for_symmetric, eslocal clearCOO_1_keep_COO_0 );
	eslocal  LoadMatrixBinInCOO(string filename, char matrix_type_G_for_general_S_for_symmetric);
	eslocal  LoadMatrixBinInCSR(string filename, char matrix_type_G_for_general_S_for_symmetric);

	void ConvertToCSR();
	void ConvertToCSR( eslocal clearCOO_1_keep_COO_0 );
	void ConvertToCSRwithSort( eslocal clearCOO_1_keep_COO_0 );

	void ConvertToCOO( eslocal clearCSR_1_keep_CSR_0 );

	void ConvertCSRToDense( eslocal clearCSR_1_keep_CSR_0 );
	void ConvertDenseToCSR( eslocal clearDense_1_keep_Dense_0 );
	void ConvertDenseToDenseFloat( eslocal clear_DoubleDense_1_keep_DoubleDense_0 );

	void DenseMatVec(SEQ_VECTOR <double> & x_in, SEQ_VECTOR <double> & y_out);
	void DenseMatVec(SEQ_VECTOR <double> & x_in, SEQ_VECTOR <double> & y_out, char T_for_transpose_N_for_not_transpose);
	void DenseMatVec(SEQ_VECTOR <double> & x_in, SEQ_VECTOR <double> & y_out, char T_for_transpose_N_for_not_transpose, eslocal x_in_vector_start_index);
	void DenseMatVec(SEQ_VECTOR <double> & x_in, SEQ_VECTOR <double> & y_out, char T_for_transpose_N_for_not_transpose, eslocal x_in_vector_start_index, eslocal y_out_vector_start_index);
	void DenseMatVec(SEQ_VECTOR <double> & x_in, SEQ_VECTOR <double> & y_out, char T_for_transpose_N_for_not_transpose, eslocal x_in_vector_start_index, eslocal y_out_vector_start_index, double beta);
	void DenseMatVec(SEQ_VECTOR <double> & x_in, SEQ_VECTOR <double> & y_out, char T_for_transpose_N_for_not_transpose, eslocal x_in_vector_start_index, eslocal y_out_vector_start_index, double beta, double alpha);

	void DenseMatMat(SparseMatrix & A_in, char trans_A, SparseMatrix & B_in, char trans_B);


	void DenseMatVecCUDA_w_Copy (SEQ_VECTOR <double> & x_in, SEQ_VECTOR <double> & y_out, char T_for_transpose_N_for_not_transpose, eslocal x_in_vector_start_index);
	void DenseMatVecCUDA_wo_Copy(SEQ_VECTOR <double> & x_in, SEQ_VECTOR <double> & y_out, char T_for_transpose_N_for_not_transpose, eslocal x_in_vector_start_index);

	void DenseMatVecCUDA_wo_Copy_start( double * x_in, double * y_out, char T_for_transpose_N_for_not_transpose, eslocal x_in_vector_start_index );
//	void DenseMatVecCUDA_shared_wo_Copy_start(double * x_in, double * y_out, char T_for_transpose_N_for_not_transpose, eslocal x_in_vector_start_index, eslocal extren_rows, eslocal extern_lda, eslocal dense_val_offset, char U_for_upper_L_for_lower);
	void DenseMatVecCUDA_shared_wo_Copy_start(double * x_in, double * y_out, char T_for_transpose_N_for_not_transpose, eslocal x_in_vector_start_index);
	void DenseMatVecCUDA_wo_Copy_sync ( );

	void DenseMatVecCUDA_wo_Copy_start_fl( float * x_in, float * y_out, char T_for_transpose_N_for_not_transpose, eslocal x_in_vector_start_index );

	void RemoveLowerDense( );

	eslocal  CopyToCUDA_Dev ();
	eslocal  CopyToCUDA_Dev_fl ();

	eslocal MallocOnCUDA_Dev (  );
	eslocal MallocOnCUDA_Dev_fl (  );
	eslocal MallocVecsOnCUDA_Dev ( );
	eslocal MallocVecsOnCUDA_Dev_fl ( );

	void CopyFromCUDA_Dev();

	void FreeFromCUDA_Dev();
	void FreeFromCUDA_Dev_fl();
	void FreeVecsFromCUDA_Dev();
	void FreeVecsFromCUDA_Dev_fl();

	void MatVecCOO(SEQ_VECTOR <double> & x_in, SEQ_VECTOR <double> & y_out);
	void MatVecCOO(SEQ_VECTOR <double> & x_in, SEQ_VECTOR <double> & y_out, char T_for_transpose_N_for_non_transpose);
	void MatVecCOO(SEQ_VECTOR <double> & x_in, SEQ_VECTOR <double> & y_out, char T_for_transpose_N_for_non_transpose, double beta);
	void MatVecCOO(SEQ_VECTOR <double> & x_in, SEQ_VECTOR <double> & y_out, char T_for_transpose_N_for_non_transpose, double beta, double alpha);

	void MatVec(SEQ_VECTOR <double> & x_in, SEQ_VECTOR <double> & y_out, char T_for_transpose_N_for_non_transpose );
	void MatVec(SEQ_VECTOR <double> & x_in, SEQ_VECTOR <double> & y_out, char T_for_transpose_N_for_non_transpose, eslocal x_in_vector_start_index, eslocal y_out_vector_start_index);
	void MatVec(SEQ_VECTOR <double> & x_in, SEQ_VECTOR <double> & y_out, char T_for_transpose_N_for_non_transpose, eslocal x_in_vector_start_index, eslocal y_out_vector_start_index, double beta);

	//void MatVecSum(vector < vector <double> > & x_in, vector <double> & y_out, char T_for_transpose_N_for_non_transpose );

	void MatMat(SparseMatrix & A_in, char MatA_T_for_transpose_N_for_non_transpose, SparseMatrix & B_in );
	void MatMatSorted(SparseMatrix & A_in, char MatA_T_for_transpose_N_for_non_transpose, SparseMatrix & B_in);

	void MatMatT(SparseMatrix & A_in, SparseMatrix & B_in);

	void MatAdd(SparseMatrix & A_in, SparseMatrix & B_in, char MatB_T_for_transpose_N_for_non_transpose, double beta);
	void MatAddInPlace(SparseMatrix & B_in, char MatB_T_for_transpose_N_for_non_transpose, double beta);

	void PrintMatSize( string Matname );


	double MatCondNumb(SparseMatrix & A_in, const std::string &str0, eslocal plot_n_first_n_last_eigenvalues,double *maxEig, int nMax_);
	void spmv_(SparseMatrix & A_in, double *x, double *Ax);
	void printMatCSR( char *str0);
	void printMatCSR2( char *str0);
	void getNullPivots(SEQ_VECTOR <eslocal> & null_pivots);
	void tridiagFromCSR( SparseMatrix & A_in, char *str0);
	double dot_e(double *x, double *y, eslocal n);

	double getNorm_K_R(SparseMatrix & K, SparseMatrix &R_in_dense_format , char str0);
	void GramSchmidtOrtho();

	void Mat_MP_Inverse(SparseMatrix &R_in, SparseMatrix &A_in);

	void MatScale(double alpha);

	void MatTranspose();
	void MatTranspose(double beta);
	void MatTranspose(SparseMatrix & A_out);
	void MatTranspose(SparseMatrix & A_out, double beta);

	void MatTransposeCOO();

	double GetMeanOfDiagonalOfSymmetricMatrix();
	double GetMaxOfDiagonalOfSymmetricMatrix();
	void   SetDiagonalOfSymmetricMatrix( double val );
	void   getSubDiagBlockmatrix( SparseMatrix & A_in, SparseMatrix & A_out, eslocal i_start, eslocal size);
	void   getSubBlockmatrix_rs( SparseMatrix & A_in, SparseMatrix & A_out, eslocal i_start, eslocal i_size,
                                          eslocal j_start, eslocal j_size);

	// Return diagonal of CSR matrix. Columns indices have to be sorted!!
	std::vector<double> getDiagonal() const;
	double getDiagonalMaximum() const;
	double getDiagonalAbsMaximum() const;

	void MatAppend(SparseMatrix & A);
	void RemoveLower();

	void CreateMatFromRowsFromMatrix         (SparseMatrix & A_in, SEQ_VECTOR <eslocal> & rows_to_add);
	void CreateMatFromRowsFromMatrix_NewSize (SparseMatrix & A_in, SEQ_VECTOR <eslocal> & rows_to_add);

	void Clear();

	eslocal MatCompare  ( SparseMatrix & A);
	eslocal MatCompareCOO(SparseMatrix & A);

	void CreateEye  ( eslocal size);
	void CreateEye  ( eslocal size, double value, eslocal offset_row, eslocal offset_col);

	void TestEye    ( eslocal size);
	void TestMatRow ( eslocal size, eslocal row_index);

	void sortInCOO();
	std::string SpyText();

//	void get_kernel_from_K();
	void get_kernel_from_K(SparseMatrix &K, SparseMatrix &regMat, SparseMatrix &KplusR,
        double &norm_KR, eslocal &defect, eslocal d_sub, size_t scSize);
  void get_kernels_from_nonsym_K(SparseMatrix &K, SparseMatrix &regMat, SparseMatrix &KplusR,
        SparseMatrix &KplusR2,
        double &norm_KR,eslocal &defect,eslocal d_sub, size_t scSize);

private:

  // get_kernel_from_K: default parameters
  //
  bool DIAGONALSCALING                                  = true;
  eslocal PERMUTVECTORACTIVE                            = 1;
  bool USE_NULL_PIVOTS_OR_S_SET                         = true;
  bool DIAGONALREGULARIZATION                           = true;
  eslocal GET_N_FIRST_AND_N_LAST_EIGENVALS_FROM_DENSE_K = 0;
  eslocal GET_N_FIRST_AND_N_LAST_EIGENVALS_FROM_DENSE_S = 0;
  eslocal PLOT_N_FIRST_N_LAST_EIGENVALUES               = 0;
  eslocal FIXING_NODES_OR_DOF                           = 0;
  eslocal DOFPERNODE                                    = 3;
  double COND_NUMB_FOR_SINGULAR_MATRIX                  = 1e13;
  eslocal CHECK_NONSING                                 = 0;
  eslocal MAX_SIZE_OF_DENSE_MATRIX_TO_GET_EIGS          = 2500;
  eslocal SC_SIZE                                       = 200;
  eslocal TWENTY                                        = 20;
  double JUMP_IN_EIGENVALUES_ALERTING_SINGULARITY       = 1.0e-5;

};

}
