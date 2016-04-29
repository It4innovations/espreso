#include "DenseSolverCUDA.h"
#include <cuda_runtime_api.h>

#define CHECK_SO(x) do {\
	cusolverStatus_t st = (x);\
	if (st != CUSOLVER_STATUS_SUCCESS) {\
		printf("API error (cusolverStatus) failed %s: %d Returned: %d\n", __FILE__, __LINE__, st);\
		switch (st) {\
			case CUSOLVER_STATUS_NOT_INITIALIZED:\
				printf("\tStatus: CUSOLVER_STATUS_NOT_INITIALIZED\n");\
				break;\
			case CUSOLVER_STATUS_ALLOC_FAILED:\
				printf("\tStatus: CUSOLVER_STATUS_ALLOC_FAILED\n");\
				break;\
			case CUSOLVER_STATUS_INVALID_VALUE:\
				printf("\tStatus: CUSOLVER_STATUS_INVALID_VALUE\n");\
				break;\
			case CUSOLVER_STATUS_ARCH_MISMATCH:\
				printf("\tStatus: CUSOLVER_STATUS_ARCH_MISMATCH\n");\
				break;\
			case CUSOLVER_STATUS_INTERNAL_ERROR:\
				printf("\tStatus: CUSOLVER_STATUS_INTERNAL_ERROR\n");\
				break;\
			case CUSOLVER_STATUS_MATRIX_TYPE_NOT_SUPPORTED:\
				printf("\tStatus: CUSOLVER_STATUS_MATRIX_TYPE_NOT_SUPPORTED\n");\
			default:\
				printf("\tUnknown status\n");\
		}\
		exit(st);\
	}\
} while(0);

#define CHECK_SO_FACT(x) do {\
	cusolverStatus_t st = (x);\
	if (st != CUSOLVER_STATUS_SUCCESS) {\
		printf("API error (cusolverStatus) failed %s: %d Returned: %d\n", __FILE__, __LINE__, st);\
		switch (st) {\
			case CUSOLVER_STATUS_NOT_INITIALIZED:\
				printf("\tStatus: CUSOLVER_STATUS_NOT_INITIALIZED\n");\
				break;\
			case CUSOLVER_STATUS_ALLOC_FAILED:\
				printf("\tStatus: CUSOLVER_STATUS_ALLOC_FAILED\n");\
				break;\
			case CUSOLVER_STATUS_INVALID_VALUE:\
				printf("\tStatus: CUSOLVER_STATUS_INVALID_VALUE\n");\
				break;\
			case CUSOLVER_STATUS_ARCH_MISMATCH:\
				printf("\tStatus: CUSOLVER_STATUS_ARCH_MISMATCH\n");\
				break;\
			case CUSOLVER_STATUS_INTERNAL_ERROR:\
				printf("\tStatus: CUSOLVER_STATUS_INTERNAL_ERROR\n");\
				break;\
			case CUSOLVER_STATUS_MATRIX_TYPE_NOT_SUPPORTED:\
				printf("\tStatus: CUSOLVER_STATUS_MATRIX_TYPE_NOT_SUPPORTED\n");\
			default:\
				printf("\tUnknown status\n");\
		}\
		return st;\
	}\
} while(0);

#define CHECK_ERR(x) do {\
	cudaError_t err = (x);\
	if (err != cudaSuccess) {\
		printf("API error (cudaError) failed %s: %d Returned: %d %s\n", __FILE__, __LINE__, err, cudaGetErrorString(err));\
		exit(err);\
	}\
} while(0);

#define CHECK_ERR_FACT(x) do {\
	cudaError_t err = (x);\
	if (err != cudaSuccess) {\
		printf("API error (cudaError) failed %s: %d Returned: %d %s\n", __FILE__, __LINE__, err, cudaGetErrorString(err));\
		return err;\
	}\
} while(0);

using namespace espreso;

DenseSolverCUDA::DenseSolverCUDA(){

	#ifdef DEBUG
		int rtVersion;
		int driverVersion;

		CHECK_ERR(cudaRuntimeGetVersion(&rtVersion));
		CHECK_ERR(cudaDriverGetVersion(&driverVersion));

		printf("CUDA runtime version: %d\nCUDA driver version: %d\n", rtVersion, driverVersion);
	#endif

// 	keep_factors=true;
// 	initialized = false;
	USE_FLOAT = false;
// 	keep_buffer = true;
	import_with_copy = false;

	cuStream 			= NULL;
	soDnHandle 			= NULL;

	D_devInfo 			= NULL;
	D_dense_values 		= NULL;
	D_dense_values_fl 	= NULL;
	D_B_dense_values 	= NULL;
	D_B_dense_values_fl = NULL;

	// if(sizeof(eslocal) == 8){
	// 	printf("64-bit integer in use with cuSolver - results may be inaccurate!\n");
	// }

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
	// m_factorized = 0;

	// Initialize cuSolver context and CUDA stream
	//CHECK_ERR(cudaSetDevice(1)); // uncomment for Espreso-WS
	while(1)
		CHECK_SO(cusolverDnCreate(&soDnHandle));
	CHECK_ERR(cudaStreamCreate(&cuStream));
	CHECK_SO(cusolverDnSetStream(soDnHandle, cuStream));
}

DenseSolverCUDA::~DenseSolverCUDA() {

	this->Clear();

	// Destroy cuSolver context and CUDA stream
	CHECK_SO(cusolverDnDestroy(soDnHandle));
	CHECK_ERR(cudaStreamDestroy(cuStream));
}


void DenseSolverCUDA::Clear() {

	MKL_INT nRhs = 1;

	CHECK_ERR(cudaFree(D_devInfo));

	if(m_dense_values_size > 0) CHECK_ERR(cudaFree(D_dense_values));
	if(m_dense_values_fl_size > 0) {
		CHECK_ERR(cudaFree(D_dense_values_fl));
		delete [] m_dense_values_fl;
	}
	if(D_B_dense_values != NULL) CHECK_ERR(cudaFree(D_B_dense_values));
	if(D_B_dense_values_fl != NULL) CHECK_ERR(cudaFree(D_B_dense_values_fl));

	m_dense_values_size = 0;
	m_dense_values_fl_size = 0;

	D_dense_values 		= NULL;
	D_dense_values_fl 	= NULL;
	m_dense_values_fl 	= NULL;
	D_B_dense_values 	= NULL;
	D_B_dense_values_fl = NULL;
}

void DenseSolverCUDA::ImportMatrix(SparseMatrix & A) {

	USE_FLOAT = false;

	m_rows	= A.rows;
	m_cols	= A.cols;
	m_nnz	= A.nnz;
	m_lda 	= A.rows;
	m_ldb 	= m_rows;

	m_dense_values_size = A.dense_values.size();
	m_dense_values_fl_size = 0;

	// cudaMalloc
	CHECK_ERR(cudaMalloc((void**)&D_devInfo, sizeof(eslocal) * 1));

	CHECK_ERR(cudaMalloc((void**)&D_dense_values, sizeof(double) * m_lda * m_rows));

	// Copy A HtoD
	CHECK_ERR(cudaMemcpy(D_dense_values, &A.dense_values[0], sizeof(double) * m_lda * m_rows, cudaMemcpyHostToDevice));

	import_with_copy = true;
}

void DenseSolverCUDA::ImportMatrix_fl(SparseMatrix & A) {

	USE_FLOAT = true;

	m_rows		= A.rows;
	m_cols		= A.cols;
	m_nnz		= A.nnz;
	m_lda		= A.rows;
	m_ldb 		= m_rows;

	m_dense_values_size = 0;
	m_dense_values_fl_size = A.dense_values.size();

	// Cast double to float
	m_dense_values_fl = new float[m_dense_values_fl_size];

	for (eslocal i = 0; i < m_dense_values_fl_size; i++)
		m_dense_values_fl[i] = (float) A.dense_values[i];

	// cudaMalloc
	CHECK_ERR(cudaMalloc((void**)&D_devInfo, sizeof(eslocal) * 1));

	CHECK_ERR(cudaMalloc((void**)&D_dense_values_fl, sizeof(float) * m_lda * m_rows));

	// Copy A HtoD
	CHECK_ERR(cudaMemcpy(D_dense_values_fl, &m_dense_values_fl[0], sizeof(float) * m_lda * m_rows, cudaMemcpyHostToDevice));

	import_with_copy = true;
}

void DenseSolverCUDA::ImportMatrix_wo_Copy(SparseMatrix & A) {
	// printf("ImportMatrix_wo_Copy: cuSolver used - Matrix will be imported to GPU!\n");
	ImportMatrix(A);
}

void DenseSolverCUDA::SetThreaded() {

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

int DenseSolverCUDA::Factorization(const std::string &str) {

	eslocal Lwork = 0;

	if (USE_FLOAT) {
		CHECK_SO_FACT(cusolverDnSpotrf_bufferSize(soDnHandle, CUBLAS_FILL_MODE_UPPER, m_rows, D_dense_values_fl, m_lda, &Lwork));

		float * D_workspace;
		CHECK_ERR_FACT(cudaMalloc((void**)&D_workspace, Lwork));

		CHECK_SO_FACT(cusolverDnSpotrf(soDnHandle, CUBLAS_FILL_MODE_UPPER, m_rows, D_dense_values_fl, m_lda, D_workspace, Lwork, D_devInfo));

		CHECK_ERR_FACT(cudaFree(D_workspace));
	} else {
		CHECK_SO_FACT(cusolverDnDpotrf_bufferSize(soDnHandle, CUBLAS_FILL_MODE_UPPER, m_rows, D_dense_values, m_lda, &Lwork));

		// Bunch-Kaufman factorization buffer
		// CHECK_SO_FACT(cusolverDnDsytrf_bufferSize(soDnHandle, m_rows, D_dense_values, m_lda, &Lwork));

		double * D_workspace;
		CHECK_ERR_FACT(cudaMalloc((void**)&D_workspace, Lwork));

		CHECK_SO_FACT(cusolverDnDpotrf(soDnHandle, CUBLAS_FILL_MODE_UPPER, m_rows, D_dense_values, m_lda, D_workspace, Lwork, D_devInfo));

		// Bunch-Kaufman factorization
		// int * D_ipiv;
		// CHECK_ERR_FACT(cudaMalloc((void**)&D_ipiv, m_rows * sizeof(int)));
		// CHECK_SO_FACT(cusolverDnDsytrf(soDnHandle, CUBLAS_FILL_MODE_UPPER, m_rows, D_dense_values, m_lda, D_ipiv, D_workspace, Lwork, D_devInfo));
		// CHECK_ERR_FACT(cudaFree(D_ipiv));		

		CHECK_ERR_FACT(cudaFree(D_workspace));
	}

	// devInfo handling
	eslocal devInfo = 0;
	CHECK_ERR_FACT(cudaMemcpy(&devInfo , D_devInfo, sizeof(eslocal) * 1, cudaMemcpyDeviceToHost));

	if (USE_FLOAT) {
		if (D_B_dense_values_fl == NULL)
			CHECK_ERR_FACT(cudaMalloc((void**)&D_B_dense_values_fl, sizeof(float) * m_ldb * m_nRhs));
	} else {
		if (D_B_dense_values == NULL)
			CHECK_ERR_FACT(cudaMalloc((void**)&D_B_dense_values, sizeof(double) * m_ldb * m_nRhs));
	}

	return devInfo;
}

void DenseSolverCUDA::Solve( SEQ_VECTOR <double> & rhs_sol) {

	m_nRhs = 1;

	if (USE_FLOAT) {
		tmp_sol_fl.resize(rhs_sol.size());

		for (eslocal i = 0; i < rhs_sol.size(); i++)
			tmp_sol_fl[i] = (float)rhs_sol[i];

		CHECK_ERR(cudaMalloc((void**)&D_B_dense_values_fl, sizeof(float) * m_ldb * m_nRhs));

		// Copy B HtoD
		CHECK_ERR(cudaMemcpy(D_B_dense_values_fl, &tmp_sol_fl[0], sizeof(float) * m_ldb * m_nRhs, cudaMemcpyHostToDevice));

		CHECK_SO(cusolverDnSpotrs(soDnHandle, CUBLAS_FILL_MODE_UPPER, m_rows, m_nRhs, D_dense_values_fl, m_lda, D_B_dense_values_fl, m_ldb, D_devInfo));

		// Copy B DtoH
		CHECK_ERR(cudaMemcpy(&tmp_sol_fl[0] , D_B_dense_values_fl, sizeof(float) * m_ldb * m_nRhs, cudaMemcpyDeviceToHost));

		for (eslocal i = 0; i < rhs_sol.size(); i++)
			rhs_sol[i] = (double)tmp_sol_fl[i];		

	} else {
		CHECK_ERR(cudaMalloc((void**)&D_B_dense_values, sizeof(double) * m_ldb * m_nRhs));

		// Copy B HtoD
		CHECK_ERR(cudaMemcpy(D_B_dense_values, &rhs_sol[0], sizeof(double) * m_ldb * m_nRhs, cudaMemcpyHostToDevice));

		CHECK_SO(cusolverDnDpotrs(soDnHandle, CUBLAS_FILL_MODE_UPPER, m_rows, m_nRhs, D_dense_values, m_lda, D_B_dense_values, m_ldb, D_devInfo));

		// Copy B DtoH
		CHECK_ERR(cudaMemcpy(&rhs_sol[0] , D_B_dense_values, sizeof(double) * m_ldb * m_nRhs, cudaMemcpyDeviceToHost));
	}
}

void DenseSolverCUDA::Solve( SEQ_VECTOR <double> & rhs, SEQ_VECTOR <double> & sol, MKL_INT n_rhs) {

	if (USE_FLOAT) {
		sol.resize(rhs.size());
		tmp_sol_fl.resize(rhs.size());

		for (eslocal i = 0; i < rhs.size(); i++)
			tmp_sol_fl[i] = (float)rhs[i];

		if (n_rhs > m_nRhs) {
			CHECK_ERR(cudaFree(D_B_dense_values_fl));
			CHECK_ERR(cudaMalloc((void**)&D_B_dense_values_fl, sizeof(float) * m_ldb * n_rhs));
		}

		// Copy B HtoD
		CHECK_ERR(cudaMemcpy(D_B_dense_values_fl, &tmp_sol_fl[0], sizeof(float) * m_ldb * n_rhs, cudaMemcpyHostToDevice));

		CHECK_SO(cusolverDnSpotrs(soDnHandle, CUBLAS_FILL_MODE_UPPER, m_rows, n_rhs, D_dense_values_fl, m_lda, D_B_dense_values_fl, m_ldb, D_devInfo));

		// Copy B DtoH
		CHECK_ERR(cudaMemcpy(&tmp_sol_fl[0] , D_B_dense_values_fl, sizeof(float) * m_ldb * n_rhs, cudaMemcpyDeviceToHost));

		for (eslocal i = 0; i < rhs.size(); i++)
			sol[i] = (double)tmp_sol_fl[i];		

	} else {
		if (n_rhs > m_nRhs) {
			CHECK_ERR(cudaFree(D_B_dense_values));
			CHECK_ERR(cudaMalloc((void**)&D_B_dense_values, sizeof(float) * m_ldb * n_rhs));
		}

		// Copy B HtoD
		CHECK_ERR(cudaMemcpy(D_B_dense_values, &rhs[0], sizeof(double) * m_ldb * n_rhs, cudaMemcpyHostToDevice));

		CHECK_SO(cusolverDnDpotrs(soDnHandle, CUBLAS_FILL_MODE_UPPER, m_rows, n_rhs, D_dense_values, m_lda, D_B_dense_values, m_ldb, D_devInfo));

		// Copy B DtoH
		CHECK_ERR(cudaMemcpy(&sol[0] , D_B_dense_values, sizeof(double) * m_ldb * n_rhs, cudaMemcpyDeviceToHost));
	}
}

void DenseSolverCUDA::Solve( SEQ_VECTOR <double> & rhs, SEQ_VECTOR <double> & sol, MKL_INT rhs_start_index, MKL_INT sol_start_index) {

	printf("DenseSolverCUDA::Solve( SEQ_VECTOR <double> & rhs, SEQ_VECTOR <double> & sol, MKL_INT rhs_start_index, MKL_INT sol_start_index) not implemented yet.\n");
	exit(1);
}
