#ifdef WIN32	 
#include "stdafx.h"
#endif

#include "utils.h"


#ifdef CUDA
	#include <cuda_runtime.h>
	#include <cublas_v2.h> 
#endif

using std::vector;
using std::cout;
using std::string; 
using std::ifstream;
using std::stringstream; 
using std::endl; 
using std::getline; 
using std::fill;


#pragma once



class SparseMatrix
{

public:
	//Constructors 
	SparseMatrix(char matrix_type_G_for_general_S_for_symmetric, string filename); 	
	SparseMatrix(char matrix_type_G_for_general_S_for_symmetric); 
	SparseMatrix();
	SparseMatrix( const SparseMatrix &A_in); 

	//Destructor 
	~SparseMatrix();

	//assignment operator 
	SparseMatrix& operator= (const SparseMatrix &A_in);

	// Variables 
	int  rows;		// number of rows 
	int  cols;		// number of columns 
	int  nnz;		// number of non zero elements  
	char type;		// 'G' for general or 'S' for symmetric 

	// Sparse COO data 
	SEQ_VECTOR <int>	I_row_indices;
	SEQ_VECTOR <int>	J_col_indices;
	SEQ_VECTOR <double> V_values; 

	// Sparse CSR data 
	SEQ_VECTOR <int>	CSR_I_row_indices;
	SEQ_VECTOR <int>	CSR_J_col_indices;
	SEQ_VECTOR <double> CSR_V_values; 

	// Dense data 
	SEQ_VECTOR <double> dense_values;
	SEQ_VECTOR <float> dense_values_fl;

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
#endif	


	//Members 

	//void LoadMatrixInCOO(string filename, char matrix_type_G_for_general_S_for_symmetric); 
	//void LoadMatrix(string filename, char matrix_type_G_for_general_S_for_symmetric); 
	int  LoadMatrixBinInCOO(string filename, char matrix_type_G_for_general_S_for_symmetric);
	int  LoadMatrixBin(string filename, char matrix_type_G_for_general_S_for_symmetric);
	int  LoadMatrixBin(string filename, char matrix_type_G_for_general_S_for_symmetric, int clearCOO_1_keep_COO_0 );

	void ConvertToCSR();
	void ConvertToCSR( int clearCOO_1_keep_COO_0 );
	void ConvertToCSRwithSort( int clearCOO_1_keep_COO_0 ); 

	void ConvertToCOO( int clearCSR_1_keep_CSR_0 );
	


	void ConvertCSRToDense( int clearCSR_1_keep_CSR_0 );
	void ConvertDenseToCSR( int clearDense_1_keep_Dense_0 );
	void ConvertDenseToDenseFloat( int clear_DoubleDense_1_keep_DoubleDense_0 );

	void DenseMatVec(SEQ_VECTOR <double> & x_in, SEQ_VECTOR <double> & y_out);
	void DenseMatVec(SEQ_VECTOR <double> & x_in, SEQ_VECTOR <double> & y_out, char T_for_transpose_N_for_not_transpose);
	void DenseMatVec(SEQ_VECTOR <double> & x_in, SEQ_VECTOR <double> & y_out, char T_for_transpose_N_for_not_transpose, int x_in_vector_start_index);

	void DenseMatVecCUDA_w_Copy (SEQ_VECTOR <double> & x_in, SEQ_VECTOR <double> & y_out, char T_for_transpose_N_for_not_transpose, int x_in_vector_start_index);
	void DenseMatVecCUDA_wo_Copy(SEQ_VECTOR <double> & x_in, SEQ_VECTOR <double> & y_out, char T_for_transpose_N_for_not_transpose, int x_in_vector_start_index);

	void DenseMatVecCUDA_wo_Copy_start( double * x_in, double * y_out, char T_for_transpose_N_for_not_transpose, int x_in_vector_start_index );
	void DenseMatVecCUDA_wo_Copy_sync ( );
	
	void DenseMatVecCUDA_wo_Copy_start_fl( float * x_in, float * y_out, char T_for_transpose_N_for_not_transpose, int x_in_vector_start_index );

	void RemoveLowerDense( ); 

	void CopyToCUDA_Dev (); 
	void CopyToCUDA_Dev_fl (); 
	void CopyFromCUDA_Dev();
	void FreeFromCUDA_Dev();

	void MatVecCOO(SEQ_VECTOR <double> & x_in, SEQ_VECTOR <double> & y_out);
	void MatVecCOO(SEQ_VECTOR <double> & x_in, SEQ_VECTOR <double> & y_out, char T_for_transpose_N_for_non_transpose);
	void MatVecCOO(SEQ_VECTOR <double> & x_in, SEQ_VECTOR <double> & y_out, char T_for_transpose_N_for_non_transpose, double beta);
	
	void MatVec(SEQ_VECTOR <double> & x_in, SEQ_VECTOR <double> & y_out, char T_for_transpose_N_for_non_transpose ); 
	void MatVec(SEQ_VECTOR <double> & x_in, SEQ_VECTOR <double> & y_out, char T_for_transpose_N_for_non_transpose, int x_in_vector_start_index, int y_out_vector_start_index);
	void MatVec(SEQ_VECTOR <double> & x_in, SEQ_VECTOR <double> & y_out, char T_for_transpose_N_for_non_transpose, int x_in_vector_start_index, int y_out_vector_start_index, double beta); 

	//void MatVecSum(vector < vector <double> > & x_in, vector <double> & y_out, char T_for_transpose_N_for_non_transpose );

	void MatMat(SparseMatrix & A_in, char MatA_T_for_transpose_N_for_non_transpose, SparseMatrix & B_in ); 
	void MatMatSorted(SparseMatrix & A_in, char MatA_T_for_transpose_N_for_non_transpose, SparseMatrix & B_in);

	void MatMatT(SparseMatrix & A_in, SparseMatrix & B_in);

	void MatAdd(SparseMatrix & A_in, SparseMatrix & B_in, char MatB_T_for_transpose_N_for_non_transpose, double beta);
	void MatAddInPlace(SparseMatrix & B_in, char MatB_T_for_transpose_N_for_non_transpose, double beta); 

	void MatScale(double alpha);

	//void MatTranspose(); 
	//void MatTranspose(double beta); 
	//void MatTranspose(SparseMatrix & A_out); 
	//void MatTranspose(SparseMatrix & A_out, double beta); 

	void MatTranspose();
	void MatTranspose(double beta);
	void MatTranspose(SparseMatrix & A_out);
	void MatTranspose(SparseMatrix & A_out, double beta);

	void MatTransposeCOO(); 

	double GetMeanOfDiagonalOfSymmetricMatrix();
	double GetMaxOfDiagonalOfSymmetricMatrix();
	
	void MatAppend(SparseMatrix & A); 
	void RemoveLower(); 

	void CreateMatFromRowsFromMatrix(SparseMatrix & A_in, SEQ_VECTOR <int> & rows_to_add);

	void Clear();

	int MatCompare  ( SparseMatrix & A);
	int MatCompareCOO(SparseMatrix & A);

	void CreateEye  ( int size); 
	void TestEye    ( int size); 
	void TestMatRow ( int size, int row_index);
	
	void sortInCOO();
};

void sortMatrixInCOO(SparseMatrix & Matrix);
void SpyText (SparseMatrix & A);
static void q_sort(SparseMatrix & Matrix, int lo, int hi );
static void q_sort_in(vector <int>    & I_row_indices, 
				   	  vector <int>    & J_col_indices, 
					  vector <double> & V_values,  
					  int lo, int hi );
