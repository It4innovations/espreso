#include "SparseSolverCUDA.h"
#include <cuda_runtime_api.h>

#define CHECK_SO(x) do {\
	cusolverStatus_t st = (x);\
	if (st != CUSOLVER_STATUS_SUCCESS) {\
		printf("API error (cusolverStatus) failed %s: %d Returned: %d\n", __FILE__, __LINE__, st);\
		printf("STATUS %d %d %d %d %d %d %d\n ", st == CUSOLVER_STATUS_SUCCESS, st == CUSOLVER_STATUS_NOT_INITIALIZED, st == CUSOLVER_STATUS_ALLOC_FAILED,\
			st == CUSOLVER_STATUS_INVALID_VALUE, st == CUSOLVER_STATUS_ARCH_MISMATCH, st == CUSOLVER_STATUS_INTERNAL_ERROR, st == CUSOLVER_STATUS_MATRIX_TYPE_NOT_SUPPORTED);\
		exit(1);\
	}\
} while(0);

#define CHECK_SP(x) do {\
	cusparseStatus_t st = (x);\
	if (st != CUSPARSE_STATUS_SUCCESS) {\
		printf("API error (cusparseStatus) failed %s: %d Returned: %d\n", __FILE__, __LINE__, st);\
		exit(1);\
	}\
} while(0);

#define CHECK_ERR(x) do {\
	cudaError_t err = (x);\
	if (err != cudaSuccess) {\
		printf("API error (cudaError) failed %s: %d Returned: %d %s\n", __FILE__, __LINE__, err, cudaGetErrorString(err));\
		exit(1);\
	}\
} while(0);

using namespace espreso;

SparseSolverCUDA::SparseSolverCUDA(){
	// printf("---SparseSolverCUDA\n");

	#ifdef DEBUG
		msglvl = 1;				/* Print statistical information in file */

		int rtVersion;
		int driverVersion;

		CHECK_ERR(cudaRuntimeGetVersion(&rtVersion));
		CHECK_ERR(cudaDriverGetVersion(&driverVersion));

		printf("CUDA runtime version: %d\nCUDA driver version: %d\n", rtVersion, driverVersion);
	#else
		msglvl = 0;
	#endif

	keep_factors=true;
	initialized = false;
	USE_FLOAT = false;
	keep_buffer = true;

	CSR_I_row_indices_size = 0;
	CSR_J_col_indices_size = 0;
	CSR_V_values_size = 0;
	CSR_V_values_fl_size = 0;
	permutation_size = 0;
	cuStream = NULL;
	soHandle = NULL;
	soInfo = NULL;
	matDescr = NULL;

	D_buffer = NULL;
	D_CSR_I_row_indices = NULL;
	D_CSR_J_col_indices = NULL;
	D_CSR_V_values = NULL;
	D_CSR_V_values_fl = NULL;
	D_rhs_sol = NULL;
	D_rhs_sol_fl = NULL;

	reorder = 2; // 0 - no reordering, 1 - symrcm, 2 - symamd

	if(sizeof(eslocal) == 8){
		printf("64-bit integer in use with cuSolver - results may be inaccurate!\n");
	}

//	/* Numbers of processors, value of OMP_NUM_THREADS */
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
	m_factorized = 0;

	// Initialize cuSolver context and CUDA stream
	//CHECK_ERR(cudaSetDevice(1)); // uncomment for Espreso-WS
	CHECK_SO(cusolverSpCreate(&soHandle));
	CHECK_ERR(cudaStreamCreate(&cuStream));
	CHECK_SO(cusolverSpSetStream(soHandle, cuStream));
}

SparseSolverCUDA::~SparseSolverCUDA() {

	this->Clear();

	// Destroy cuSolver context and CUDA stream
	CHECK_SO(cusolverSpDestroy(soHandle));
	CHECK_ERR(cudaStreamDestroy(cuStream));
}


void SparseSolverCUDA::Clear() {
	// printf("---Clear\n");

	if (USE_FLOAT){
		if(D_rhs_sol_fl != NULL){
			CHECK_ERR(cudaFree(D_rhs_sol_fl));
			D_rhs_sol_fl = NULL;
		}
	} else {
		if (D_rhs_sol != NULL){
			CHECK_ERR(cudaFree(D_rhs_sol));
			D_rhs_sol = NULL;
		}
	}

	if (initialized == true) {
		MKL_INT nRhs = 1;

		// Destroys factorization and matrix settings
		if(soInfo != NULL){
			CHECK_SO(cusolverSpDestroyCsrcholInfo(soInfo));
			soInfo = NULL;
		}

		// if (D_buffer != NULL){
		// 	CHECK_ERR(cudaFree(D_buffer));
		// 	D_buffer = NULL;
		// }

		initialized = false;
	}
	// else if (keep_buffer == true && D_buffer != NULL) {
	// 	CHECK_ERR(cudaFree(D_buffer));
	// 	D_buffer = NULL;
	// }
	if (D_buffer != NULL){
		CHECK_ERR(cudaFree(D_buffer));
		D_buffer = NULL;
	}

	if (matDescr != NULL){
		CHECK_SP(cusparseDestroyMatDescr(matDescr));
		matDescr = NULL;
	}

	if (CSR_I_row_indices_size > 0){
		// delete [] CSR_I_row_indices; // TODO probably delete
		// CSR_I_row_indices = NULL; // TODO probably delete

		CHECK_ERR(cudaFree(D_CSR_I_row_indices));
		D_CSR_I_row_indices = NULL;

		CSR_I_row_indices_size = 0;
	}
	if (CSR_J_col_indices_size > 0)	{
		// delete [] CSR_J_col_indices; // TODO probably delete
		// CSR_J_col_indices = NULL; // TODO probably delete

		CHECK_ERR(cudaFree(D_CSR_J_col_indices));
		D_CSR_J_col_indices = NULL;

		CSR_J_col_indices_size = 0;
	}
	if (CSR_V_values_size > 0){
		// delete [] CSR_V_values;	// TODO probably delete
		// CSR_V_values = NULL;	// TODO probably delete

		CHECK_ERR(cudaFree(D_CSR_V_values));
		D_CSR_V_values = NULL;

		CSR_V_values_size = 0;
	}
	if (CSR_V_values_fl_size > 0){
		CHECK_ERR(cudaFree(D_CSR_V_values_fl));
		D_CSR_V_values_fl = NULL;

		CSR_V_values_fl_size   = 0;
	}
	if (permutation_size > 0){
		delete [] permutation;
		permutation = NULL;

		permutation_size = 0;
	}
}

void SparseSolverCUDA::ReorderMatrix(SparseMatrix & A) {
		eslocal i;

		MKL_INT rows_r 						= A.rows;
		MKL_INT nnz_r					 	= A.nnz;
		MKL_INT CSR_I_row_indices_size_r	= A.CSR_I_row_indices.size();
		MKL_INT CSR_J_col_indices_size_r	= A.CSR_J_col_indices.size();
		MKL_INT CSR_V_values_size_r			= A.CSR_V_values.size();

		// printf("rows: %d nnz: %d i: %d j: %d v: %d \n", rows_r, nnz_r, CSR_I_row_indices_size_r, CSR_J_col_indices_size_r, A.CSR_V_values.size());
		// SpyText(A);

		permutation_size = rows_r;
		permutation = new int[permutation_size];

		SEQ_VECTOR<double> CSR_V_values_reordered(CSR_V_values_size_r);

		#if INT_WIDTH == 64
			SEQ_VECTOR<int> tmp_CSR_I_row_indices(A.CSR_I_row_indices.begin(), A.CSR_I_row_indices.end());
			SEQ_VECTOR<int> tmp_CSR_J_col_indices(A.CSR_J_col_indices.begin(), A.CSR_J_col_indices.end());
		#endif

		// Create permutation vector
		if(reorder == 1){
			// printf("Symrcm reordering method performed.\n");
			#if INT_WIDTH == 64
				CHECK_SO(cusolverSpXcsrsymrcmHost(soHandle, rows_r, nnz_r, matDescr,  &tmp_CSR_I_row_indices.front(),
				 &tmp_CSR_J_col_indices.front(), permutation));
			#else
				CHECK_SO(cusolverSpXcsrsymrcmHost(soHandle, rows_r, nnz_r, matDescr,  &A.CSR_I_row_indices.front(),
				 &A.CSR_J_col_indices.front(), permutation));
			#endif
		}
		else if(reorder == 2){
			// printf("Symamd reordering method performed.\n");
			#if INT_WIDTH == 64
				CHECK_SO(cusolverSpXcsrsymamdHost(soHandle, rows_r, nnz_r, matDescr, &tmp_CSR_I_row_indices.front(),
				 &tmp_CSR_J_col_indices.front(), permutation));
			#else
				CHECK_SO(cusolverSpXcsrsymamdHost(soHandle, rows_r, nnz_r, matDescr, &A.CSR_I_row_indices.front(),
				 &A.CSR_J_col_indices.front(), permutation));
			#endif
		}

		#ifdef DEBUG
			printf("Permutation vector: ");
			for (i = 0; i < rows_r; ++i)
				printf("%d ", permutation[i]);
			printf("\n");
		#endif

		// Allocate a buffer for reordering
		size_t bufferSizeInBytes = 0;
		#if INT_WIDTH == 64
			CHECK_SO(cusolverSpXcsrperm_bufferSizeHost(soHandle, rows_r, rows_r, nnz_r, matDescr,
		 	 &tmp_CSR_I_row_indices.front(), &tmp_CSR_J_col_indices.front(), permutation, permutation, &bufferSizeInBytes));
		#else
			CHECK_SO(cusolverSpXcsrperm_bufferSizeHost(soHandle, rows_r, rows_r, nnz_r, matDescr,
		 	 &A.CSR_I_row_indices.front(), &A.CSR_J_col_indices.front(), permutation, permutation, &bufferSizeInBytes));
		#endif

		#ifdef DEBUG
			printf("---Permutation bufferSizeInBytes: %zu B\n", bufferSizeInBytes);
		#endif

		void * buffer = malloc(bufferSizeInBytes);

		int * map = new int[nnz_r];
		for (i = 0; i < nnz_r; ++i){
			map[i] = i;
		}

		// Perform the matrix reordering (get a map)
		#if INT_WIDTH == 64
			CHECK_SO(cusolverSpXcsrpermHost(soHandle, rows_r, rows_r, nnz_r, matDescr,
			 &tmp_CSR_I_row_indices.front(), &tmp_CSR_J_col_indices.front(), permutation, permutation, map, buffer));

			for (i = 0; i < CSR_I_row_indices_size_r; ++i)
				A.CSR_I_row_indices[i] = tmp_CSR_I_row_indices[i];

			for (i = 0; i < CSR_J_col_indices_size_r; ++i)
				A.CSR_J_col_indices[i] = tmp_CSR_J_col_indices[i];
		#else
			CHECK_SO(cusolverSpXcsrpermHost(soHandle, rows_r, rows_r, nnz_r, matDescr,
			 &A.CSR_I_row_indices.front(), &A.CSR_J_col_indices.front(), permutation, permutation, map, buffer));
		#endif

		// Get reordered values
		for (i = 0; i < nnz_r; ++i){
			CSR_V_values_reordered[i] = A.CSR_V_values[map[i]];
		}

		A.CSR_V_values = CSR_V_values_reordered;

		delete [] map;
		free(buffer);
}

void SparseSolverCUDA::ImportMatrix(SparseMatrix & A_in) {
	// printf("---ImportMatrix\n");

	USE_FLOAT = false;

	CHECK_SP(cusparseCreateMatDescr(&matDescr));
	cusparseSetMatType(matDescr, CUSPARSE_MATRIX_TYPE_GENERAL);
	cusparseSetMatIndexBase(matDescr, CUSPARSE_INDEX_BASE_ONE);

	SparseMatrix* A;
	SparseMatrix A_sym;

	eslocal i;

	#ifdef DEBUG
		#if INT_WIDTH == 64
			printf("Input matrix: rows: %ld nnz: %ld i: %zu j: %zu v: %zu \n", A_in.rows, A_in.nnz, A_in.CSR_I_row_indices.size(), A_in.CSR_J_col_indices.size(), A_in.CSR_V_values.size());
		#else
			printf("Input matrix: rows: %d nnz: %d i: %zu j: %zu v: %zu \n", A_in.rows, A_in.nnz, A_in.CSR_I_row_indices.size(), A_in.CSR_J_col_indices.size(), A_in.CSR_V_values.size());
		#endif

		A_in.SpyText();
	#endif

	if(reorder > 0)
	{
		A_sym = A_in;
		// Extend upper part to full matrix for reordering
		if(A_in.type == 'S'){
			A_sym.type = 'G';
			A_sym.SetDiagonalOfSymmetricMatrix(0.0);
			A_sym.MatTranspose();
			A_sym.MatAddInPlace(A_in,'N',1.0);
		}

		ReorderMatrix(A_sym);

		#ifdef DEBUG
			A_sym.SpyText();
		#endif

		A_sym.RemoveLower();
		A_sym.type = 'S';

		A = &A_sym;

		#ifdef DEBUG
			(*A).SpyText();
		#endif
	}
	else{
		// Import only lower part if stored as G
		if(A_in.type == 'G'){
			A_sym = A_in;
			A_sym.RemoveLower();
			A_sym.type = 'S';
			A = &A_sym;

			printf("Matrix transformed from type G to type S internally.\n");
		} else {
			A = &A_in;
		}

		// printf("No reordering performed.\n");
	}

	A->MatTranspose();

	#ifdef DEBUG
		#if INT_WIDTH == 64
			printf("Transposed matrix: rows: %ld nnz: %ld i: %zu j: %zu v: %zu \n", A->rows, A->nnz, A->CSR_I_row_indices.size(), A->CSR_J_col_indices.size(), A->CSR_V_values.size());
		#else
			printf("Transposed matrix: rows: %d nnz: %d i: %zu j: %zu v: %zu \n", A->rows, A->nnz, A->CSR_I_row_indices.size(), A->CSR_J_col_indices.size(), A->CSR_V_values.size());
		#endif

		(*A).SpyText();
	#endif

	rows	= A->rows;
	cols	= A->cols;
	nnz		= A->nnz;
	m_Kplus_size = A->rows;

	CSR_I_row_indices_size = A->CSR_I_row_indices.size();
	CSR_J_col_indices_size = A->CSR_J_col_indices.size();
	CSR_V_values_size	   = A->CSR_V_values.size();

	// ORIG
	// CSR_I_row_indices = new int[CSR_I_row_indices_size];
	// CSR_J_col_indices = new int[CSR_J_col_indices_size];
	// CSR_V_values	  = new double  [CSR_V_values_size];

	CHECK_ERR(cudaMalloc ((void**)&D_CSR_I_row_indices, sizeof(int)*CSR_I_row_indices_size));
	CHECK_ERR(cudaMalloc ((void**)&D_CSR_J_col_indices, sizeof(int)*CSR_J_col_indices_size));
	CHECK_ERR(cudaMalloc ((void**)&D_CSR_V_values, sizeof(double)*CSR_V_values_size));

	CHECK_ERR(cudaMemcpy(D_CSR_V_values, &A->CSR_V_values.front(), sizeof(double)*CSR_V_values_size, cudaMemcpyHostToDevice));

	#if INT_WIDTH == 64
		SEQ_VECTOR<int> tmp2_CSR_I_row_indices(A->CSR_I_row_indices.begin(), A->CSR_I_row_indices.end());
		SEQ_VECTOR<int> tmp2_CSR_J_col_indices(A->CSR_J_col_indices.begin(), A->CSR_J_col_indices.end());

		CHECK_ERR(cudaMemcpy(D_CSR_I_row_indices, &tmp2_CSR_I_row_indices.front(), sizeof(int)*CSR_I_row_indices_size, cudaMemcpyHostToDevice));
		CHECK_ERR(cudaMemcpy(D_CSR_J_col_indices, &tmp2_CSR_J_col_indices.front(), sizeof(int)*CSR_J_col_indices_size, cudaMemcpyHostToDevice));
	#else
		CHECK_ERR(cudaMemcpy(D_CSR_I_row_indices, &A->CSR_I_row_indices.front(), sizeof(int)*CSR_I_row_indices_size, cudaMemcpyHostToDevice));
		CHECK_ERR(cudaMemcpy(D_CSR_J_col_indices, &A->CSR_J_col_indices.front(), sizeof(int)*CSR_J_col_indices_size, cudaMemcpyHostToDevice));
	#endif

	A->MatTranspose();

	// ORIG
	// copy(A_in.CSR_I_row_indices.begin(), A_in.CSR_I_row_indices.end(), CSR_I_row_indices);
	// copy(A_in.CSR_J_col_indices.begin(), A_in.CSR_J_col_indices.end(), CSR_J_col_indices);
	// copy(A_in.CSR_V_values     .begin(), A_in.CSR_V_values     .end(), CSR_V_values);

	import_with_copy = true;
}

void SparseSolverCUDA::ImportMatrix_fl(SparseMatrix & A_in) {
	// std::printf("---ImportMatrix_fl\n");

	USE_FLOAT = true;

	CHECK_SP(cusparseCreateMatDescr(&matDescr));
	cusparseSetMatType(matDescr, CUSPARSE_MATRIX_TYPE_GENERAL);
	cusparseSetMatIndexBase(matDescr, CUSPARSE_INDEX_BASE_ONE);

	SparseMatrix* A;
	SparseMatrix A_sym;

	eslocal i;

	#ifdef DEBUG
		#if INT_WIDTH == 64
			printf("Input matrix: rows: %ld nnz: %ld i: %zu j: %zu v: %zu \n", A_in.rows, A_in.nnz, A_in.CSR_I_row_indices.size(), A_in.CSR_J_col_indices.size(), A_in.CSR_V_values.size());
		#else
			printf("Input matrix: rows: %d nnz: %d i: %zu j: %zu v: %zu \n", A_in.rows, A_in.nnz, A_in.CSR_I_row_indices.size(), A_in.CSR_J_col_indices.size(), A_in.CSR_V_values.size());
		#endif

		A_in.SpyText();
	#endif

	if(reorder > 0)
	{
		A_sym = A_in;
		// Extend upper part to full matrix for reordering
		if(A_in.type == 'S'){
			A_sym.type = 'G';
			A_sym.SetDiagonalOfSymmetricMatrix(0.0);
			A_sym.MatTranspose();
			A_sym.MatAddInPlace(A_in,'N',1.0);
		}

		ReorderMatrix(A_sym);

		#ifdef DEBUG
			A_sym.SpyText();
		#endif
		A_sym.RemoveLower();
		A_sym.type = 'S';

		A = &A_sym;

		#ifdef DEBUG
			(*A).SpyText();
		#endif
	} else {
		// Import only lower part if stored as G
		if(A_in.type == 'G'){
			A_sym = A_in;
			A_sym.RemoveLower();
			A_sym.type = 'S';
			A = &A_sym;

			printf("Matrix transformed from type G to type S internally.\n");
		} else {
			A = &A_in;
		}
	}

	A->MatTranspose();

	#ifdef DEBUG
		#if INT_WIDTH == 64
			printf("Transposed matrix: rows: %ld nnz: %ld i: %zu j: %zu v: %zu \n", A->rows, A->nnz, A->CSR_I_row_indices.size(), A->CSR_J_col_indices.size(), A->CSR_V_values.size());
		#else
			printf("Transposed matrix: rows: %d nnz: %d i: %zu j: %zu v: %zu \n", A->rows, A->nnz, A->CSR_I_row_indices.size(), A->CSR_J_col_indices.size(), A->CSR_V_values.size());
		#endif

		(*A).SpyText();
	#endif

	rows	= A->rows;
	cols	= A->cols;
	nnz		= A->nnz;
	m_Kplus_size = A->rows;

	CSR_I_row_indices_size = A->CSR_I_row_indices.size();
	CSR_J_col_indices_size = A->CSR_J_col_indices.size();
	CSR_V_values_fl_size   = A->CSR_V_values.size();

	// ORIG
	// CSR_I_row_indices = new int[CSR_I_row_indices_size];
	// CSR_J_col_indices = new int[CSR_J_col_indices_size];
	// CSR_V_values	  = new double  [CSR_V_values_size];

	// Cast double to float
	SEQ_VECTOR<float> CSR_V_values_fl(A->CSR_V_values.begin(), A->CSR_V_values.end());

	CHECK_ERR(cudaMalloc ((void**)&D_CSR_I_row_indices, sizeof(int)*CSR_I_row_indices_size));
	CHECK_ERR(cudaMalloc ((void**)&D_CSR_J_col_indices, sizeof(int)*CSR_J_col_indices_size));
	CHECK_ERR(cudaMalloc ((void**)&D_CSR_V_values_fl, sizeof(float)*CSR_V_values_fl_size));

	CHECK_ERR(cudaMemcpy(D_CSR_V_values_fl, &CSR_V_values_fl.front(), sizeof(float)*CSR_V_values_fl_size, cudaMemcpyHostToDevice));

	#if INT_WIDTH == 64
		SEQ_VECTOR<int> tmp_CSR_I_row_indices(A->CSR_I_row_indices.begin(), A->CSR_I_row_indices.end());
		SEQ_VECTOR<int> tmp_CSR_J_col_indices(A->CSR_J_col_indices.begin(), A->CSR_J_col_indices.end());

		CHECK_ERR(cudaMemcpy(D_CSR_I_row_indices, &tmp_CSR_I_row_indices.front(), sizeof(int)*CSR_I_row_indices_size, cudaMemcpyHostToDevice));
		CHECK_ERR(cudaMemcpy(D_CSR_J_col_indices, &tmp_CSR_J_col_indices.front(), sizeof(int)*CSR_J_col_indices_size, cudaMemcpyHostToDevice));
	#else
		CHECK_ERR(cudaMemcpy(D_CSR_I_row_indices, &A->CSR_I_row_indices.front(), sizeof(int)*CSR_I_row_indices_size, cudaMemcpyHostToDevice));
		CHECK_ERR(cudaMemcpy(D_CSR_J_col_indices, &A->CSR_J_col_indices.front(), sizeof(int)*CSR_J_col_indices_size, cudaMemcpyHostToDevice));
	#endif

	A->MatTranspose();

	// ORIG
	// copy(A_in.CSR_I_row_indices.begin(), A_in.CSR_I_row_indices.end(), CSR_I_row_indices);
	// copy(A_in.CSR_J_col_indices.begin(), A_in.CSR_J_col_indices.end(), CSR_J_col_indices);
	// copy(A_in.CSR_V_values     .begin(), A_in.CSR_V_values     .end(), CSR_V_values);

	import_with_copy = true;
}

void SparseSolverCUDA::ImportMatrix_wo_Copy(SparseMatrix & A) {
	// printf("ImportMatrix_wo_Copy: cuSolver used - Matrix will be imported to GPU!\n");
	ImportMatrix(A);
}

void SparseSolverCUDA::SetThreaded() {

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

int SparseSolverCUDA::Factorization(const std::string &str) {
	// printf("---Factorization\n");

	// Keeps factorization
	CHECK_SO(cusolverSpCreateCsrcholInfo(&soInfo));

	if(matDescr == NULL){
		printf("Factorization error: cuSolver matrix description structure doesn't exist!\n");
		exit(1);
	}

	CHECK_SO(cusolverSpXcsrcholAnalysis(soHandle, rows, nnz, matDescr, D_CSR_I_row_indices, D_CSR_J_col_indices, soInfo));

	internalDataInBytes = 0;
	workspaceInBytes = 0;

	if (USE_FLOAT) {
		CHECK_SO(cusolverSpScsrcholBufferInfo(soHandle, rows, nnz, matDescr, D_CSR_V_values_fl, D_CSR_I_row_indices, D_CSR_J_col_indices, soInfo, &internalDataInBytes, &workspaceInBytes));

		if (D_buffer == NULL)
			CHECK_ERR(cudaMalloc ((void**)&D_buffer, workspaceInBytes));

		CHECK_SO(cusolverSpScsrcholFactor(soHandle, rows, nnz, matDescr, D_CSR_V_values_fl, D_CSR_I_row_indices, D_CSR_J_col_indices, soInfo, D_buffer));

		if (!keep_buffer){
			CHECK_ERR(cudaFree(D_buffer));
			D_buffer = NULL;
		}
	} else {
		CHECK_SO(cusolverSpDcsrcholBufferInfo(soHandle, rows, nnz, matDescr, D_CSR_V_values, D_CSR_I_row_indices, D_CSR_J_col_indices, soInfo, &internalDataInBytes, &workspaceInBytes));

		if (D_buffer == NULL)
			CHECK_ERR(cudaMalloc ((void**)&D_buffer, workspaceInBytes));

		CHECK_SO(cusolverSpDcsrcholFactor(soHandle, rows, nnz, matDescr, D_CSR_V_values, D_CSR_I_row_indices, D_CSR_J_col_indices, soInfo, D_buffer));

		if (!keep_buffer){
			CHECK_ERR(cudaFree(D_buffer));
			D_buffer = NULL;
		}
	}

	#ifdef DEBUG
		printf("---internalDataInBytes: %zu B\n", internalDataInBytes);
		printf("---workspaceInBytes: %zu B\n", workspaceInBytes);

		printf ("\nFactorization completed ... ");
		printf("I: %d, J: %d, V: %d, RHS: %d, celkem: %d\n", (rows+1)*sizeof(int), nnz*sizeof(int), nnz*sizeof(float), rows*sizeof(float),
			(rows+1)*sizeof(int)+ nnz*sizeof(int)+ nnz*sizeof(float)+ rows*sizeof(float));
	#endif

	m_factorized = 1; // used only in constructor and factorization
	initialized = true;

	if (USE_FLOAT) {
		D_rhs_sol_fl = NULL;
		CHECK_ERR(cudaMalloc ((void**)&D_rhs_sol_fl, sizeof(float)*m_Kplus_size));
		rhs_sol_fl.resize(m_Kplus_size);
	} else {
		D_rhs_sol = NULL;
		CHECK_ERR(cudaMalloc ((void**)&D_rhs_sol, sizeof(double)*m_Kplus_size));
		if(reorder)
			rhs_sol_reordered.resize(m_Kplus_size);
	}

	return 0;
}

void SparseSolverCUDA::Solve( SEQ_VECTOR <double> & rhs_sol) {
	// printf("---SOLVE1 %d \n", USE_FLOAT);

	if (!initialized){
		std::stringstream ss;
		ss << "Solve -> rank: " << config::MPIrank;
		Factorization(ss.str());
	}

	int n_rhs    = 1;
	eslocal i;

	/* -------------------------------------------------------------------- */
	/* .. Back substitution and iterative refinement. */
	/* -------------------------------------------------------------------- */
	if (USE_FLOAT){
		// Reordering RHS
		if(reorder){
			for (i = 0; i < permutation_size; ++i){
				rhs_sol_fl[i] = (float)rhs_sol[permutation[i]];
			}
		} else {
			for (i = 0; i < rows; ++i){
				rhs_sol_fl[i] = (float)rhs_sol[i];
			}
		}

		CHECK_ERR(cudaMemcpy(D_rhs_sol_fl, &rhs_sol_fl.front(), sizeof(float)*rows, cudaMemcpyHostToDevice));

		if(D_buffer == NULL)
			CHECK_ERR(cudaMalloc ((void**)&D_buffer, workspaceInBytes));

		CHECK_SO(cusolverSpScsrcholSolve(soHandle, rows, D_rhs_sol_fl, D_rhs_sol_fl, soInfo, D_buffer));

		if (!keep_buffer){
			CHECK_ERR(cudaFree(D_buffer));
			D_buffer = NULL;
		}

		CHECK_ERR(cudaMemcpy(&rhs_sol_fl.front(), D_rhs_sol_fl, sizeof(float)*rows, cudaMemcpyDeviceToHost));

		// Reorder and cast back the solution
		if(reorder){
			for (i = 0; i < permutation_size; ++i){
	        	rhs_sol[permutation[i]] = (double)rhs_sol_fl[i];
	        }
		} else {
			for (i = 0; i < rows; i++){
				rhs_sol[i] = (double)rhs_sol_fl[i];
			}
		}

	} else {
		// Reordering RHS
		if(reorder){
			// rhs_sol_reordered.resize(permutation_size);

			for (i = 0; i < permutation_size; ++i){
				rhs_sol_reordered[i] = rhs_sol[permutation[i]];
			}
		}
		CHECK_ERR(cudaMemcpy(D_rhs_sol, reorder ? &rhs_sol_reordered.front() : &rhs_sol.front(), sizeof(double)*rows, cudaMemcpyHostToDevice));

		if(D_buffer == NULL)
			CHECK_ERR(cudaMalloc ((void**)&D_buffer, workspaceInBytes));

		CHECK_SO(cusolverSpDcsrcholSolve(soHandle, rows, D_rhs_sol, D_rhs_sol, soInfo, D_buffer));

		if (!keep_buffer){
			CHECK_ERR(cudaFree(D_buffer));
			D_buffer = NULL;
		}

		CHECK_ERR(cudaMemcpy(reorder ? &rhs_sol_reordered.front() : &rhs_sol.front(), D_rhs_sol, sizeof(double)*rows, cudaMemcpyDeviceToHost));

        // Reorder back the solution
        for (i = 0; i < permutation_size; ++i){
        	rhs_sol[permutation[i]] = rhs_sol_reordered[i];
        }
	}

	if (!keep_factors) {
		/* -------------------------------------------------------------------- */
		/* .. Termination and release of memory. */
		/* -------------------------------------------------------------------- */
		int nRhs = 1;
		initialized = false;
		if (MPIrank == 0) printf(".");

		if(soInfo != NULL){
			CHECK_SO(cusolverSpDestroyCsrcholInfo(soInfo));
			soInfo = NULL;
		}

		if(USE_FLOAT){
			CHECK_ERR(cudaFree(D_rhs_sol_fl));
			D_rhs_sol_fl = NULL;
		}else{
			CHECK_ERR(cudaFree(D_rhs_sol));
			D_rhs_sol = NULL;
		}

		// printf("Factors freed\n");
	}

	if (!keep_buffer){
		CHECK_ERR(cudaFree(D_buffer));
		D_buffer = NULL;
	}
}

void SparseSolverCUDA::Solve( SEQ_VECTOR <double> & rhs, SEQ_VECTOR <double> & sol, MKL_INT n_rhs) {
	// std::printf("---SOLVE2 n_rhs=%d \n", n_rhs);

	if (!initialized){
		std::stringstream ss;
		ss << "Solve -> rank: " << config::MPIrank;
		Factorization(ss.str());
	}

	eslocal i;

	// for (int i = 0; i < n_rhs; ++i)
	// {
	// 	SEQ_VECTOR <double> one_rhs_sol(rows);
	// 	copy(rhs.begin() + i * rows, rhs.begin() + rows + i * rows, one_rhs_sol.begin());
	// 	Solve(one_rhs_sol);
	// 	copy(one_rhs_sol.begin(), one_rhs_sol.end(), sol.begin() + i * rows);
	// }

	eslocal total_size = rows * n_rhs; // if reorder -> rows == permutation_size

	if (USE_FLOAT){
		if(m_Kplus_size < total_size){
			rhs_sol_fl.resize(total_size);
			CHECK_ERR(cudaFree(D_rhs_sol_fl));
			D_rhs_sol_fl = NULL;
			CHECK_ERR(cudaMalloc ((void**)&D_rhs_sol_fl, sizeof(float) * total_size));
		}

		// Reordering RHS
		if(reorder){
			for (i = 0; i < total_size; ++i){
				rhs_sol_fl[i] = (float)rhs[permutation[i % permutation_size] + (i / permutation_size) * permutation_size];
			}
		} else {
			for (i = 0; i < total_size; ++i){
				rhs_sol_fl[i] = (float)rhs[i];
			}
		}

		CHECK_ERR(cudaMemcpy(D_rhs_sol_fl, &rhs_sol_fl.front(), sizeof(float) * total_size, cudaMemcpyHostToDevice));

		if(D_buffer == NULL)
			CHECK_ERR(cudaMalloc ((void**)&D_buffer, workspaceInBytes));

		for (i = 0; i < n_rhs; ++i){
			CHECK_SO(cusolverSpScsrcholSolve(soHandle, rows, D_rhs_sol_fl + i * rows, D_rhs_sol_fl + i * rows, soInfo, D_buffer));
		}

		if (!keep_buffer){
			CHECK_ERR(cudaFree(D_buffer));
			D_buffer = NULL;
		}

		CHECK_ERR(cudaMemcpy(&rhs_sol_fl.front(), D_rhs_sol_fl, sizeof(float) * total_size, cudaMemcpyDeviceToHost));

		// Reorder and cast back the solution
		if(reorder){
			for (i = 0; i < total_size; ++i){
	        	sol[permutation[i % permutation_size] + (i / permutation_size) * permutation_size] = (double)rhs_sol_fl[i];
	        }
		} else {
			for (i = 0; i < total_size; i++){
				sol[i] = (double)rhs_sol_fl[i];
			}
		}
	} else {
		if(m_Kplus_size < total_size){
			CHECK_ERR(cudaFree(D_rhs_sol));
			D_rhs_sol = NULL;
			CHECK_ERR(cudaMalloc ((void**)&D_rhs_sol, sizeof(double) * total_size));
		}

		// Reordering RHS
		if(reorder){
			if(m_Kplus_size < total_size)
				rhs_sol_reordered.resize(total_size);

			for (i = 0; i < total_size; ++i){
				rhs_sol_reordered[i] = rhs[ permutation[i % permutation_size] + (i/permutation_size)*permutation_size];
			}
		}

		CHECK_ERR(cudaMemcpy(D_rhs_sol, reorder ? &rhs_sol_reordered.front() : &rhs.front(), sizeof(double) * total_size, cudaMemcpyHostToDevice));

		if(D_buffer == NULL)
			CHECK_ERR(cudaMalloc ((void**)&D_buffer, workspaceInBytes));

		for (i = 0; i < n_rhs; ++i){
			CHECK_SO(cusolverSpDcsrcholSolve(soHandle, rows, D_rhs_sol + i * rows, D_rhs_sol + i * rows, soInfo, D_buffer));
		}

		if (!keep_buffer){
			CHECK_ERR(cudaFree(D_buffer));
			D_buffer = NULL;
		}

		CHECK_ERR(cudaMemcpy(reorder ? &rhs_sol_reordered.front() : &sol.front(), D_rhs_sol, sizeof(double) * total_size, cudaMemcpyDeviceToHost));

		if(reorder){
			// Reorder back the solution
	        for (i = 0; i < total_size; ++i){
	        	sol[permutation[i % permutation_size] + (i / permutation_size) * permutation_size] = rhs_sol_reordered[i];
	        }
		}
	}

	/* -------------------------------------------------------------------- */
	/* .. Back substitution and iterative refinement. */
	/* -------------------------------------------------------------------- */

	if (!keep_factors) {
		/* -------------------------------------------------------------------- */
		/* .. Termination and release of memory. */
		/* -------------------------------------------------------------------- */
		MKL_INT nRhs = 1;
		initialized = false;
		if (MPIrank == 0) printf(".");

		if(soInfo != NULL){
			CHECK_SO(cusolverSpDestroyCsrcholInfo(soInfo));
			soInfo = NULL;
		}

		if(USE_FLOAT){
			CHECK_ERR(cudaFree(D_rhs_sol_fl));
			D_rhs_sol_fl = NULL;
		}else{
			CHECK_ERR(cudaFree(D_rhs_sol));
			D_rhs_sol = NULL;
		}

		// printf("Factors freed\n");
	}

	if (!keep_buffer){
		CHECK_ERR(cudaFree(D_buffer));
		D_buffer = NULL;
	}
}

void SparseSolverCUDA::Solve( SEQ_VECTOR <double> & rhs, SEQ_VECTOR <double> & sol, MKL_INT rhs_start_index, MKL_INT sol_start_index) {
	// std::printf("---SOLVE3 n_rhs=%d \n", m_nRhs);

	if (!initialized){
		std::stringstream ss;
		ss << "Solve -> rank: " << config::MPIrank;
		Factorization(ss.str());
	}

	eslocal i;

	// SEQ_VECTOR <double> one_rhs_sol(rows);
	// copy(&rhs[rhs_start_index], &rhs[rhs_start_index + rows], one_rhs_sol.begin());
	// Solve(one_rhs_sol);
	// copy(one_rhs_sol.begin(), one_rhs_sol.end(), &sol[sol_start_index]);

	if (USE_FLOAT){
		// Reordering RHS
		if(reorder){
			for (i = 0; i < permutation_size; ++i){
				rhs_sol_fl[i] = (float)rhs[rhs_start_index + permutation[i]];
			}
		} else {
			for (i = 0; i < rows; ++i){
				rhs_sol_fl[i] = (float)rhs[rhs_start_index + i];
			}
		}

		CHECK_ERR(cudaMemcpy(D_rhs_sol_fl, &rhs_sol_fl.front(), sizeof(float)*rows, cudaMemcpyHostToDevice));

		if(D_buffer == NULL)
			CHECK_ERR(cudaMalloc ((void**)&D_buffer, workspaceInBytes));

		CHECK_SO(cusolverSpScsrcholSolve(soHandle, rows, D_rhs_sol_fl, D_rhs_sol_fl, soInfo, D_buffer));

		if (!keep_buffer){
			CHECK_ERR(cudaFree(D_buffer));
			D_buffer = NULL;
		}

		CHECK_ERR(cudaMemcpy(&rhs_sol_fl.front(), D_rhs_sol_fl, sizeof(float)*rows, cudaMemcpyDeviceToHost));

		// Reorder and cast back the solution
		if(reorder){
			for (i = 0; i < permutation_size; ++i){
	        	sol[sol_start_index + permutation[i]] = (double)rhs_sol_fl[i];
	        }
		} else {
			for (i = 0; i < rows; i++){
				sol[sol_start_index + i] = (double)rhs_sol_fl[i];
			}
		}
	} else {
		// Reordering RHS
		if(reorder){
			for (i = 0; i < permutation_size; ++i){
				rhs_sol_reordered[i] = rhs[rhs_start_index + permutation[i]];
			}
		}

		CHECK_ERR(cudaMemcpy(D_rhs_sol, reorder ? &rhs_sol_reordered.front() : &rhs[rhs_start_index], sizeof(double)*rows, cudaMemcpyHostToDevice));

		if(D_buffer == NULL)
			CHECK_ERR(cudaMalloc ((void**)&D_buffer, workspaceInBytes));

		CHECK_SO(cusolverSpDcsrcholSolve(soHandle, rows, D_rhs_sol, D_rhs_sol, soInfo, D_buffer));

		if (!keep_buffer){
			CHECK_ERR(cudaFree(D_buffer));
			D_buffer = NULL;
		}

		CHECK_ERR(cudaMemcpy(reorder ? &rhs_sol_reordered.front() : &sol[sol_start_index], D_rhs_sol, sizeof(double)*rows, cudaMemcpyDeviceToHost));

        // Reorder back the solution
        for (i = 0; i < permutation_size; ++i){
        	sol[sol_start_index + permutation[i]] = rhs_sol_reordered[i];
        }
	}

	if (!keep_factors) {
		/* -------------------------------------------------------------------- */
		/* .. Termination and release of memory. */
		/* -------------------------------------------------------------------- */
		initialized = false;
		if (MPIrank == 0) printf(".");

		if(soInfo != NULL){
			CHECK_SO(cusolverSpDestroyCsrcholInfo(soInfo));
			soInfo = NULL;
		}

		if(USE_FLOAT){
			CHECK_ERR(cudaFree(D_rhs_sol_fl));
			D_rhs_sol_fl = NULL;
		}else{
			CHECK_ERR(cudaFree(D_rhs_sol));
			D_rhs_sol = NULL;
		}

		// printf("Factors freed\n");
	}

	if (!keep_buffer){
		CHECK_ERR(cudaFree(D_buffer));
		D_buffer = NULL;
	}
}


void SparseSolverCUDA::SolveMat_Sparse( SparseMatrix & A) {
	SolveMat_Sparse(A, A);
};

void SparseSolverCUDA::SolveMat_Sparse( SparseMatrix & A_in, SparseMatrix & B_out) {
	SolveMat_Sparse(A_in, B_out, 'T');
};

void SparseSolverCUDA::SolveMat_Sparse( SparseMatrix & A_in, SparseMatrix & B_out, char T_for_input_matrix_is_transposed_N_input_matrix_is_NOT_transposed ) {
	// printf("---SolveMat_Sparse\n");

	//SolveMat_Dense(A_in, B_out);

	if (!initialized){
		std::stringstream ss;
		ss << "Solve -> rank: " << config::MPIrank;
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
		phase = -1;			/* Release internal memory. */
		MKL_INT nRhs = 1;
		double ddum;			/* Double dummy */
		MKL_INT idum;			/* Integer dummy. */
		// PARDISO (pt, &maxfct, &mnum, &mtype, &phase,
		// 		&rows, &ddum, CSR_I_row_indices, CSR_J_col_indices, &idum, &nRhs,
		// 		iparm, &msglvl, &ddum, &ddum, &error);
		initialized = false;

		if(soInfo != NULL){
			CHECK_SO(cusolverSpDestroyCsrcholInfo(soInfo));
			soInfo = NULL;
		}
	}
}


void SparseSolverCUDA::SolveMat_Dense( SparseMatrix & A ) {
	SolveMat_Dense(A, A);
}

void SparseSolverCUDA::SolveMat_Dense( SparseMatrix & A_in, SparseMatrix & B_out ) {
	// printf("---SolveMat_Dense\n");
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
void SparseSolverCUDA::SolveMatF( SparseMatrix & A_in, SparseMatrix & B_out, bool isThreaded ) {

	printf("Method SolveMatF is not implemented yet.\n");
	exit(1);

	// /* Internal solver memory pointer pt, */
	// /* 32-bit: int pt[64]; 64-bit: long int pt[64] */
	// /* or void *pt[64] should be OK on both architectures */
	// void *pt[64];

	// /* Pardiso control parameters. */
	// MKL_INT iparm[64];
	// MKL_INT maxfct, mnum, phase, error; //, msglvl;
	// /* Auxiliary variables. */
	// MKL_INT i;
	// double ddum;			/* Double dummy */
	// MKL_INT idum;			/* Integer dummy. */

	// /* -------------------------------------------------------------------- */
	// /* .. Setup Pardiso control parameters. */
	// /* -------------------------------------------------------------------- */
	// for (i = 0; i < 64; i++) {
	// 	iparm[i] = 0;
	// }

	// MKL_INT mtype = 2;

	// iparm[0] = 1;		/* No solver default */
	// iparm[1] = 2;		/* Fill-in reordering from METIS */
	// 					/* Numbers of processors, value of OMP_NUM_THREADS */
	// //iparm[2] = 8;		/* Not used in MKL PARDISO */ // TODO: zjistit co to je pro MKL to bylo 0


	// if (isThreaded) {
	// 	/* Numbers of processors, value of OMP_NUM_THREADS */
	// 	int num_procs;
	// 	char * var = getenv("SOLVER_NUM_THREADS");
	//     if(var != NULL)
	//     	sscanf( var, "%d", &num_procs );
	// 	else {
	//     	printf("Set environment SOLVER_NUM_THREADS to 1");
	//         exit(1);
	// 	}

	//     iparm[2] = num_procs;
	// } else {
	// 	iparm[2] = 1;
	// }



	// iparm[3] = 0;		/* No iterative-direct algorithm */
	// iparm[4] = 0;		/* No user fill-in reducing permutation */
	// iparm[5] = 0;		/* Write solution into x */
	// iparm[6] = 0;		/* Not in use */
	// iparm[7] = 0;		/* Max numbers of iterative refinement steps */
	// iparm[8] = 0;		/* Not in use */
	// iparm[9] = 13;		/* Perturb the pivot elements with 1E-13 */
	// iparm[10] = 0;		/* Use nonsymmetric permutation and scaling MPS */
	// iparm[11] = 0;		/* Not in use */
	// iparm[12] = 0;		/* Maximum weighted matching algorithm is switched-off */
	// 					/* (default for symmetric). Try iparm[12] = 1 in case of inappropriate accuracy */
	// iparm[13] = 0;		/* Output: Number of perturbed pivots */
	// iparm[14] = 0;		/* Not in use */
	// iparm[15] = 0;		/* Not in use */
	// iparm[16] = 0;		/* Not in use */
	// iparm[17] = -1;		/* Output: Number of nonzeros in the factor LU */
	// iparm[18] = -1;		/* Output: Mflops for LU factorization */
	// iparm[19] = 0;		/* Output: Numbers of CG Iterations */

	// maxfct = 1;			/* Maximum number of numerical factorizations. */
	// mnum   = 1;			/* Which factorization to use. */
	// //msglvl = 0;			/* Supress printing statistical information */
	// error  = 0;			/* Initialize error flag */

	// /* -------------------------------------------------------------------- */
	// /* .. Initialize the internal solver memory pointer. This is only */
	// /* necessary for the FIRST call of the PARDISO solver. */
	// /* -------------------------------------------------------------------- */

	// for (i = 0; i < 64; i++) pt[i] = 0;

	// MKL_INT job[8];

	// MKL_INT m		= A_in.rows;
	// MKL_INT n		= A_in.cols;
	// MKL_INT nRhs	= A_in.cols;
	// MKL_INT lda     = m;
	// MKL_INT info;

	// SEQ_VECTOR<double>  sol  (m * n, 0);

	// bool clear_dense = false;

	// //TODO: na konci se musi smazat !!
	// if (A_in.dense_values.size() == 0) {
	// 	A_in.ConvertCSRToDense(0);
	// 	clear_dense = true;
	// }

	// SEQ_VECTOR<MKL_INT> perm (A_in.dense_values.size() , 0);
	// for (MKL_INT ii = 0; ii < A_in.dense_values.size(); ii++)
	// 	if (A_in.dense_values[ii] != 0.0)
	// 		perm[ii] = 1;

	// phase = 13;
	// PARDISO (pt, &maxfct, &mnum, &mtype, &phase,
	// 	    &rows, CSR_V_values, CSR_I_row_indices, CSR_J_col_indices, &perm[0], &nRhs, iparm, &msglvl, &A_in.dense_values[0], &sol[0], &error);


	// if (error != 0)
	// {
	// 	printf ("\nERROR during the solution of the system : %d", error);
	// 	exit (1);
	// } else {
	// 	initialized = true;
	// }

	// MKL_INT x = 0;

	// // Convert solution matrix (SOL) to sparse format - find nnz step
	// job[0] = 0; // If job(1)=0, the rectangular matrix A is converted to the CSR format;
	// job[1] = 1; // if job(2)=1, one-based indexing for the rectangular matrix A is used.
	// job[2] = 1; // if job(3)=1, one-based indexing for the matrix in CSR format is used.
	// job[3] = 2; // If job(4)=2, adns is a whole matrix A.

	// job[4] = 1; // job(5)=nzmax - maximum number of the non-zero elements allowed if job(1)=0.
	// job[5] = 0; // job(6) - job indicator for conversion to CSR format.
	// 			// If job(6)=0, only array ia is generated for the output storage.
	// 			// If job(6)>0, arrays acsr, ia, ja are generated for the output storage.
	// job[6] = 0; //
	// job[7] = 0; //

	// B_out.CSR_I_row_indices.resize(m + 1);
	// B_out.CSR_J_col_indices.resize(1);
	// B_out.CSR_V_values.     resize(1);

	// mkl_ddnscsr (
	// 	job,
	// 	&m, &n,
	// 	&sol[0], &lda,
	// 	&B_out.CSR_V_values[0], &B_out.CSR_J_col_indices[0], &B_out.CSR_I_row_indices[0],
	// 	&info);

	// // Convert solution matrix (SOL) to sparse format - convert step
	// MKL_INT nnzmax = B_out.CSR_I_row_indices[m];//-1; POZOR

	// B_out.CSR_J_col_indices.resize(nnzmax);
	// B_out.CSR_V_values.     resize(nnzmax);

	// job[4] = nnzmax; // job(5) = nzmax - maximum number of the non-zero elements allowed if job(1)=0.
	// job[5] = 1; // job(6) - job indicator for conversion to CSR format.
	// 			// If job(6)=0, only array ia is generated for the output storage.
	// 			// If job(6)>0, arrays acsr, ia, ja are generated for the output storage.

	// mkl_ddnscsr (
	// 	job,
	// 	&m, &n,
	// 	&sol[0], &lda,
	// 	&B_out.CSR_V_values[0], &B_out.CSR_J_col_indices[0], &B_out.CSR_I_row_indices[0],
	// 	&info);


	// // Setup parameters for output matrix
	// B_out.cols	= A_in.cols;
	// B_out.rows	= A_in.rows;
	// B_out.nnz	= B_out.CSR_V_values.size();
	// B_out.type	= 'G';

	// SEQ_VECTOR<double>().swap( sol );

	// /* -------------------------------------------------------------------- */
	// /* .. Termination and release of memory. */
	// /* -------------------------------------------------------------------- */
	// phase = -1;			/* Release internal memory. */
	// PARDISO (pt, &maxfct, &mnum, &mtype, &phase,
	// 	&rows, &ddum, CSR_I_row_indices, CSR_J_col_indices, &idum, &nRhs,
	// 	iparm, &msglvl, &ddum, &ddum, &error);

	// //SEQ_VECTOR<double>().swap( A_in.dense_values );

	// initialized = false;

	// if (clear_dense) {
	// 	SEQ_VECTOR<double>().swap( A_in.dense_values );
	// }
}

void SparseSolverCUDA::Create_SC( SparseMatrix & SC_out, MKL_INT sc_size, bool isThreaded ) {

	printf("Method Create_SC is not implemented yet.\n");
	exit(1);
// 	/* Internal solver memory pointer pt, */
// 	/* 32-bit: int pt[64]; 64-bit: long int pt[64] */
// 	/* or void *pt[64] should be OK on both architectures */
// 	void *pt[64];

// 	/* Pardiso control parameters. */
// 	MKL_INT 	iparm[64];
// 	double  dparm[65];
// 	MKL_INT 	maxfct, mnum, phase, error;

// 	/* Auxiliary variables. */
// 	MKL_INT 	i;
// 	double 	ddum;			/* Double dummy */
// 	MKL_INT 	idum;			/* Integer dummy. */
// 	MKL_INT 	solver;

// 	/* -------------------------------------------------------------------- */
// 	/* .. Setup Pardiso control parameters. */
// 	/* -------------------------------------------------------------------- */
// 	for (i = 0; i < 64; i++) {
// 		iparm[i] = 0;
// 	}

// 	/* -------------------------------------------------------------------- */
// 	/* .. Initialize the internal solver memory pointer. This is only */
// 	/* necessary for the FIRST call of the PARDISO solver. */
// 	/* -------------------------------------------------------------------- */
// 	for (i = 0; i < 64; i++)
// 		pt[i] = 0;

// 	MKL_INT 	mtype = 2;

// 	/* Numbers of processors, value of OMP_NUM_THREADS */
// 	if (isThreaded) {
// 		/* Numbers of processors, value of OMP_NUM_THREADS */
// 		int num_procs;
// 		char * var = getenv("SOLVER_NUM_THREADS");
// 	    if(var != NULL)
// 	    	sscanf( var, "%d", &num_procs );
// 		else {
// 	    	printf("Set environment SOLVER_NUM_THREADS to 1");
// 	        exit(1);
// 		}

// 	    iparm[2] = num_procs;
// 	} else {
// 		iparm[2] = 1;
// 	}

// //	iparm[0] = 1;		/* No solver default */
// //	iparm[1] = 2;		/* Fill-in reordering from METIS */
// //	iparm[2]			/* Numbers of processors, value of OMP_NUM_THREADS */
// //	iparm[2] = 8;		/* Not used in MKL PARDISO */
// //	iparm[3] = 0;		/* No iterative-direct algorithm */
// //	iparm[4] = 0;		/* No user fill-in reducing permutation */
// //	iparm[5] = 0;		/* Write solution into x */
// //	iparm[6] = 0;		/* Not in use */
// //	iparm[7] = 0;		/* Max numbers of iterative refinement steps */
// //	iparm[8] = 0;		/* Not in use */
// //	iparm[9] = 13;		/* Perturb the pivot elements with 1E-13 */
// //	iparm[10] = 0;		/* Use nonsymmetric permutation and scaling MPS */
// //	iparm[11] = 0;		/* Not in use */
// //	iparm[12] = 0;		/* Maximum weighted matching algorithm is switched-off */
// //						/* (default for symmetric). Try iparm[12] = 1 in case of inappropriate accuracy */
// //	iparm[13] = 0;		/* Output: Number of perturbed pivots */
// //	iparm[14] = 0;		/* Not in use */
// //	iparm[15] = 0;		/* Not in use */
// //	iparm[16] = 0;		/* Not in use */
// //	iparm[17] = -1;		/* Output: Number of nonzeros in the factor LU */
// //	iparm[18] = -1;		/* Output: Mflops for LU factorization */
// //	iparm[19] = 0;		/* Output: Numbers of CG Iterations */
// //
// //	maxfct = 1;			/* Maximum number of numerical factorizations. */
// //	mnum   = 1;			/* Which factorization to use. */
// //	//msglvl = 0;			/* Supress printing statistical information */
// //	error  = 0;			/* Initialize error flag */




//     iparm[1-1] = 1;         /* No solver default */
//     iparm[2-1] = 2;         /* Fill-in reordering from METIS */
//     iparm[10-1] = 8; //13   /* Perturb the pivot elements with 1E-13 */
//     iparm[11-1] = 0;        /* Use nonsymmetric permutation and scaling MPS */
//     iparm[13-1] = 0;         Maximum weighted matching algorithm is switched-off (default for symmetric). Try iparm[12] = 1 in case of inappropriate accuracy
//     iparm[14-1] = 0;        /* Output: Number of perturbed pivots */
//     iparm[18-1] = -1;       /* Output: Number of nonzeros in the factor LU */
//     iparm[19-1] = -1;       /* Output: Mflops for LU factorization */
//     iparm[36-1] = 1;        /* Use Schur complement */

//     maxfct = 1;           /* Maximum number of numerical factorizations. */
//     mnum = 1;             /* Which factorization to use. */
//     //msglvl = 1;           /* Print statistical information in file */
//     error = 0;            /* Initialize error flag */

//     /* -------------------------------------------------------------------- */
//     /* .. Reordering and Symbolic Factorization. This step also allocates   */
//     /* all memory that is necessary for the factorization.                  */
//     /* -------------------------------------------------------------------- */

//     std::vector <MKL_INT> perm (rows ,0);
//     for (MKL_INT i = rows - sc_size; i < rows; i++)
//     	perm[i] = 1;

//     MKL_INT nrhs = 0;

// //    phase = 11;
// //    PARDISO (pt, &maxfct, &mnum, &mtype, &phase,
// //        		&K_sc1.rows,
// //				&K_sc1.CSR_V_values[0], &K_sc1.CSR_I_row_indices[0], &K_sc1.CSR_J_col_indices[0],
// //				&perm[0], &nrhs, iparm, &msglvl, &ddum, &ddum, &error);
// //
// //    if ( error != 0 )
// //    {
// //    	printf ("\nERROR during symbolic factorization: %d", error);
// //    	exit(1);
// //    }


//     /* -------------------------------------------------------------------- */
//     /* .. Numerical factorization. */
//     /* -------------------------------------------------------------------- */

// 	SC_out.dense_values.resize(sc_size * sc_size);

//     phase = 12;
//     PARDISO (pt, &maxfct, &mnum, &mtype, &phase,
// 			&rows,
// 			&CSR_V_values[0], &CSR_I_row_indices[0], &CSR_J_col_indices[0],
// 			&perm[0], &nrhs,
// 			iparm, &msglvl, &ddum, &SC_out.dense_values[0], &error);

//     //for (MKL_INT i = 0; i < SC_out.dense_values.size(); i++)
//     //	SC_out.dense_values[i] = (-1.0)*SC_out.dense_values[i];

//     if ( error != 0 )
// 	{
// 		printf ("\nERROR during numerical factorization: %d", error);
// 		exit (2);
// 	} else {
// 		initialized = true;
// 	}

// 	/* -------------------------------------------------------------------- */
// 	/* .. Termination and release of memory. */
// 	/* -------------------------------------------------------------------- */
// 	phase = -1;           /* Release internal memory. */
// 	PARDISO (pt, &maxfct, &mnum, &mtype, &phase,
// 			&rows, &ddum, &CSR_I_row_indices[0], &CSR_J_col_indices[0], &idum, &nrhs,
// 			 iparm, &msglvl, &ddum, &ddum, &error);

// 	initialized = false;

//     /* -------------------------------------------------------------------- */
//     /* ..  allocate memory for the Schur-complement and copy it there.      */
//     /* -------------------------------------------------------------------- */
//     MKL_INT nonzeros_S = iparm[38];

//     SC_out.cols = sc_size;
//     SC_out.rows = sc_size;
//     SC_out.type = 'G';

//     SC_out.ConvertDenseToCSR(1);
//     SC_out.RemoveLower();
//     SC_out.type = 'S';

//     if (msglvl == 1)
//     	SpyText(SC_out);

// //    if (generate_symmetric_sc_1_generate_general_sc_0 == 1) {
// //    	SC_out.RemoveLower();
// //    }

}

void SparseSolverCUDA::Create_SC_w_Mat( SparseMatrix & K_in, SparseMatrix & B_in, SparseMatrix & SC_out,
								    bool isThreaded, MKL_INT generate_symmetric_sc_1_generate_general_sc_0 ) {

	printf("Method Create_SC_w_Mat is not implemented yet.\n");
	exit(1);

// 	//int msglvl = 0;

// 	// *** Prepare matrix

// 	SparseMatrix K_sc1;
// 	SparseMatrix Sc_eye;
// 	SparseMatrix K_b_tmp;

// 	K_b_tmp = B_in;
// 	K_b_tmp.MatTranspose();

// 	Sc_eye.CreateEye(K_b_tmp.rows, 0.0, 0, K_b_tmp.cols);

// 	K_sc1 = K_in;
// 	K_sc1.MatTranspose();
// 	K_sc1.MatAppend(K_b_tmp);
// 	K_sc1.MatTranspose();
// 	K_sc1.MatAppend(Sc_eye);


// 	// *** END - Prepare matrix




// 	/* Internal solver memory pointer pt, */
// 	/* 32-bit: int pt[64]; 64-bit: long int pt[64] */
// 	/* or void *pt[64] should be OK on both architectures */
// 	void *pt[64];

// 	/* Pardiso control parameters. */
// 	MKL_INT 	iparm[64];
// 	double  	dparm[65];
// 	MKL_INT 	maxfct, mnum, phase, error;

// 	/* Auxiliary variables. */
// 	MKL_INT 	i;
// 	double 		ddum;			/* Double dummy */
// 	MKL_INT 	idum;			/* Integer dummy. */
// 	MKL_INT 	solver;

// 	/* -------------------------------------------------------------------- */
// 	/* .. Setup Pardiso control parameters. */
// 	/* -------------------------------------------------------------------- */
// 	for (i = 0; i < 64; i++) {
// 		iparm[i] = 0;
// 	}

// 	/* -------------------------------------------------------------------- */
// 	/* .. Initialize the internal solver memory pointer. This is only */
// 	/* necessary for the FIRST call of the PARDISO solver. */
// 	/* -------------------------------------------------------------------- */
// 	for (i = 0; i < 64; i++)
// 		pt[i] = 0;

// 	MKL_INT 	mtype = 2;

// 	/* Numbers of processors, value of OMP_NUM_THREADS */
// 	if (isThreaded) {
// 		/* Numbers of processors, value of OMP_NUM_THREADS */
// 		MKL_INT num_procs;
// 		char * var = getenv("SOLVER_NUM_THREADS");
// 	    if(var != NULL)
// 	    	sscanf( var, "%d", &num_procs );
// 		else {
// 	    	printf("Set environment SOLVER_NUM_THREADS to 1");
// 	        exit(1);
// 		}

// 	    iparm[2] = num_procs;
// 	} else {
// 		iparm[2] = 1;
// 	}

// //	iparm[0] = 1;		/* No solver default */
// //	iparm[1] = 2;		/* Fill-in reordering from METIS */
// //	iparm[2]			/* Numbers of processors, value of OMP_NUM_THREADS */
// //	iparm[2] = 8;		/* Not used in MKL PARDISO */
// //	iparm[3] = 0;		/* No iterative-direct algorithm */
// //	iparm[4] = 0;		/* No user fill-in reducing permutation */
// //	iparm[5] = 0;		/* Write solution into x */
// //	iparm[6] = 0;		/* Not in use */
// //	iparm[7] = 0;		/* Max numbers of iterative refinement steps */
// //	iparm[8] = 0;		/* Not in use */
// //	iparm[9] = 13;		/* Perturb the pivot elements with 1E-13 */
// //	iparm[10] = 0;		/* Use nonsymmetric permutation and scaling MPS */
// //	iparm[11] = 0;		/* Not in use */
// //	iparm[12] = 0;		/* Maximum weighted matching algorithm is switched-off */
// //						/* (default for symmetric). Try iparm[12] = 1 in case of inappropriate accuracy */
// //	iparm[13] = 0;		/* Output: Number of perturbed pivots */
// //	iparm[14] = 0;		/* Not in use */
// //	iparm[15] = 0;		/* Not in use */
// //	iparm[16] = 0;		/* Not in use */
// //	iparm[17] = -1;		/* Output: Number of nonzeros in the factor LU */
// //	iparm[18] = -1;		/* Output: Mflops for LU factorization */
// //	iparm[19] = 0;		/* Output: Numbers of CG Iterations */
// //
// //	maxfct = 1;			/* Maximum number of numerical factorizations. */
// //	mnum   = 1;			/* Which factorization to use. */
// //	//msglvl = 0;			/* Supress printing statistical information */
// //	error  = 0;			/* Initialize error flag */




//     iparm[1-1] = 1;         /* No solver default */
//     iparm[2-1] = 2;         /* Fill-in reordering from METIS */
//     iparm[10-1] = 8; //13   /* Perturb the pivot elements with 1E-13 */
//     iparm[11-1] = 0;        /* Use nonsymmetric permutation and scaling MPS */
//     iparm[13-1] = 0;         Maximum weighted matching algorithm is switched-off (default for symmetric). Try iparm[12] = 1 in case of inappropriate accuracy
//     iparm[14-1] = 0;        /* Output: Number of perturbed pivots */
//     iparm[18-1] = -1;       /* Output: Number of nonzeros in the factor LU */
//     iparm[19-1] = -1;       /* Output: Mflops for LU factorization */
//     iparm[36-1] = 1;        /* Use Schur complement */

//     maxfct = 1;           /* Maximum number of numerical factorizations. */
//     mnum = 1;             /* Which factorization to use. */
//     //msglvl = 1;           /* Print statistical information in file */
//     error = 0;            /* Initialize error flag */

//     /* -------------------------------------------------------------------- */
//     /* .. Reordering and Symbolic Factorization. This step also allocates   */
//     /* all memory that is necessary for the factorization.                  */
//     /* -------------------------------------------------------------------- */

//     std::vector <MKL_INT> perm (K_sc1.rows,0);
//     for (MKL_INT i = K_in.rows; i < K_sc1.rows; i++)
//     	perm[i] = 1;

//     MKL_INT nrhs = 0;

// //    phase = 11;
// //    PARDISO (pt, &maxfct, &mnum, &mtype, &phase,
// //        		&K_sc1.rows,
// //				&K_sc1.CSR_V_values[0], &K_sc1.CSR_I_row_indices[0], &K_sc1.CSR_J_col_indices[0],
// //				&perm[0], &nrhs, iparm, &msglvl, &ddum, &ddum, &error);
// //
// //    if ( error != 0 )
// //    {
// //    	printf ("\nERROR during symbolic factorization: %d", error);
// //    	exit(1);
// //    }


//     /* -------------------------------------------------------------------- */
//     /* .. Numerical factorization. */
//     /* -------------------------------------------------------------------- */

// 	SC_out.dense_values.resize(K_b_tmp.rows * K_b_tmp.rows);

//     phase = 12;
//     PARDISO (pt, &maxfct, &mnum, &mtype, &phase,
// 			&K_sc1.rows,
// 			&K_sc1.CSR_V_values[0], &K_sc1.CSR_I_row_indices[0], &K_sc1.CSR_J_col_indices[0],
// 			&perm[0], &nrhs,
// 			iparm, &msglvl, &ddum, &SC_out.dense_values[0], &error);

//     for (MKL_INT i = 0; i < SC_out.dense_values.size(); i++)
//     	SC_out.dense_values[i] = (-1.0)*SC_out.dense_values[i];

//     if ( error != 0 )
// 	{
// 		printf ("\nERROR during numerical factorization: %d", error);
// 		exit (2);
// 	} else {
// 		initialized = true;
// 	}

// 	/* -------------------------------------------------------------------- */
// 	/* .. Termination and release of memory. */
// 	/* -------------------------------------------------------------------- */
// 	phase = -1;           /* Release internal memory. */
// 	PARDISO (pt, &maxfct, &mnum, &mtype, &phase,
// 			&K_sc1.rows, &ddum, &K_sc1.CSR_I_row_indices[0], &K_sc1.CSR_J_col_indices[0], &idum, &nrhs,
// 			 iparm, &msglvl, &ddum, &ddum, &error);

// 	initialized = false;

//     /* -------------------------------------------------------------------- */
//     /* ..  allocate memory for the Schur-complement and copy it there.      */
//     /* -------------------------------------------------------------------- */
//     MKL_INT nonzeros_S = iparm[38];

//     SC_out.cols = K_b_tmp.rows;
//     SC_out.rows = K_b_tmp.rows;
//     SC_out.type = 'G';

// //    SC_out.ConvertDenseToCSR(1);

// //    if (msglvl == 1)
// //    	SpyText(SC_out);

//     if (generate_symmetric_sc_1_generate_general_sc_0 == 1) {
//     	//SC_out.RemoveLower();
//     	SC_out.RemoveLowerDense();
//     	SC_out.type = 'S';
//     }


// //    // Finalize shape of the SC
// //    if (generate_symmetric_sc_1_generate_general_sc_0 == 0) {
// //
// //		SC_out.type = 'G';
// //
// //		SparseMatrix SC_tmp;
// //		SC_tmp = SC_out;
// //		SC_tmp.SetDiagonalOfSymmetricMatrix(0.0);
// //		SC_tmp.MatTranspose();
// //
// //		SC_out.MatAddInPlace(SC_tmp,'N',1.0);
// //
// //    }
// //
// //	SC_out.MatScale(-1.0);
// //	//SC_out.ConvertCSRToDense(0);

}

void SparseSolverCUDA::Create_non_sym_SC_w_Mat( SparseMatrix & K_in, SparseMatrix & B1_in, SparseMatrix & B0_in, SparseMatrix & SC_out, bool isThreaded, MKL_INT generate_symmetric_sc_1_generate_general_sc_0 ) {

	printf("Method Create_non_sym_SC_w_Mat is not implemented yet.\n");
	exit(1);

// 	//int msglvl = 0;

// 	//SpyText(K_in);
// 	//SpyText(B1_in);
// 	//SpyText(B0_in);


// 	SparseMatrix K;

// 	// Create "non-symmetric" K matrix;
// 	if (K_in.type = 'S') {
// 		K = K_in;
// 		K.SetDiagonalOfSymmetricMatrix(0.0);
// 		K.MatTranspose();
// 		K.MatAddInPlace(K_in, 'N', 1.0);
// 		K.type = 'G';
// 	} else {
// 		// Pozor - not optimal, zbytecna kopie matice K
// 		K = K_in;
// 	}


// 	//SpyText(K);

// 	// *** Prepare matrix

// 	SparseMatrix K_sc1;
// 	SparseMatrix Sc_eye;
// 	SparseMatrix K_b_tmp;

// 	K_b_tmp = B1_in;
// 	K_b_tmp.MatTranspose();

// 	Sc_eye.CreateEye(K_b_tmp.rows, 0.0, 0, K_b_tmp.cols);

// 	K.MatTranspose();
// 	K.MatAppend(K_b_tmp);
// 	K.MatTranspose();

//     //SpyText(K);

// 	SparseMatrix K_b0_tmp;
// 	K_b0_tmp = B0_in;
// 	K_b0_tmp.MatTranspose();
// 	K_b0_tmp.ConvertToCOO(1);
// 	K_b0_tmp.rows = B1_in.cols;
// 	K_b0_tmp.cols = B1_in.cols + K_in.cols;
// 	K_b0_tmp.ConvertToCSRwithSort(1);

//     //SpyText(K_b0_tmp);

// 	K_b0_tmp.MatAddInPlace(Sc_eye,'N',1.0);

//     //SpyText(K_b0_tmp);

// 	K.MatAppend(K_b0_tmp);

//     K_sc1 = K;

// 	// *** END - Prepare matrix
// 	//SpyText(K_sc1);


// 	/* Internal solver memory pointer pt, */
// 	/* 32-bit: int pt[64]; 64-bit: long int pt[64] */
// 	/* or void *pt[64] should be OK on both architectures */
// 	void *pt[64];

// 	/* Pardiso control parameters. */
// 	MKL_INT 	iparm[64];
// 	double  dparm[65];
// 	MKL_INT 	maxfct, mnum, phase, error;

// 	/* Auxiliary variables. */
// 	MKL_INT 	i;
// 	double 	ddum;			/* Double dummy */
// 	MKL_INT 	idum;			/* Integer dummy. */
// 	MKL_INT 	solver;

// 	/* -------------------------------------------------------------------- */
// 	/* .. Setup Pardiso control parameters. */
// 	/* -------------------------------------------------------------------- */
// 	for (i = 0; i < 64; i++) {
// 		iparm[i] = 0;
// 	}

// 	/* -------------------------------------------------------------------- */
// 	/* .. Initialize the internal solver memory pointer. This is only */
// 	/* necessary for the FIRST call of the PARDISO solver. */
// 	/* -------------------------------------------------------------------- */
// 	for (i = 0; i < 64; i++)
// 		pt[i] = 0;

// 	MKL_INT 	mtype = 11;

// 	/* Numbers of processors, value of OMP_NUM_THREADS */
// 	if (isThreaded) {
// 		/* Numbers of processors, value of OMP_NUM_THREADS */
// 		int num_procs;
// 		char * var = getenv("SOLVER_NUM_THREADS");
// 	    if(var != NULL)
// 	    	sscanf( var, "%d", &num_procs );
// 		else {
// 	    	printf("Set environment SOLVER_NUM_THREADS to 1");
// 	        exit(1);
// 		}

// 	    iparm[2] = num_procs;
// 	} else {
// 		iparm[2] = 1;
// 	}

// //	iparm[0] = 1;		/* No solver default */
// //	iparm[1] = 2;		/* Fill-in reordering from METIS */
// //	iparm[2]			/* Numbers of processors, value of OMP_NUM_THREADS */
// //	iparm[2] = 8;		/* Not used in MKL PARDISO */
// //	iparm[3] = 0;		/* No iterative-direct algorithm */
// //	iparm[4] = 0;		/* No user fill-in reducing permutation */
// //	iparm[5] = 0;		/* Write solution into x */
// //	iparm[6] = 0;		/* Not in use */
// //	iparm[7] = 0;		/* Max numbers of iterative refinement steps */
// //	iparm[8] = 0;		/* Not in use */
// //	iparm[9] = 13;		/* Perturb the pivot elements with 1E-13 */
// //	iparm[10] = 0;		/* Use nonsymmetric permutation and scaling MPS */
// //	iparm[11] = 0;		/* Not in use */
// //	iparm[12] = 0;		/* Maximum weighted matching algorithm is switched-off */
// //						/* (default for symmetric). Try iparm[12] = 1 in case of inappropriate accuracy */
// //	iparm[13] = 0;		/* Output: Number of perturbed pivots */
// //	iparm[14] = 0;		/* Not in use */
// //	iparm[15] = 0;		/* Not in use */
// //	iparm[16] = 0;		/* Not in use */
// //	iparm[17] = -1;		/* Output: Number of nonzeros in the factor LU */
// //	iparm[18] = -1;		/* Output: Mflops for LU factorization */
// //	iparm[19] = 0;		/* Output: Numbers of CG Iterations */
// //
// //	maxfct = 1;			/* Maximum number of numerical factorizations. */
// //	mnum   = 1;			/* Which factorization to use. */
// //	//msglvl = 0;			/* Supress printing statistical information */
// //	error  = 0;			/* Initialize error flag */




//     iparm[1-1] = 1;         /* No solver default */
//     iparm[2-1] = 2;         /* Fill-in reordering from METIS */
//     iparm[10-1] = 8; //13   /* Perturb the pivot elements with 1E-13 */
//     iparm[11-1] = 0;        /* Use nonsymmetric permutation and scaling MPS */
//     iparm[13-1] = 0;        /* Maximum weighted matching algorithm is switched-off (default for symmetric). Try iparm[12] = 1 in case of inappropriate accuracy */
//     iparm[14-1] = 0;        /* Output: Number of perturbed pivots */
//     iparm[18-1] = -1;       /* Output: Number of nonzeros in the factor LU */
//     iparm[19-1] = -1;       /* Output: Mflops for LU factorization */
//     iparm[36-1] = 1;        /* Use Schur complement */

//     maxfct = 1;           /* Maximum number of numerical factorizations. */
//     mnum = 1;             /* Which factorization to use. */
//     //msglvl = 1;           /* Print statistical information in file */
//     error = 0;            /* Initialize error flag */

//     /* -------------------------------------------------------------------- */
//     /* .. Reordering and Symbolic Factorization. This step also allocates   */
//     /* all memory that is necessary for the factorization.                  */
//     /* -------------------------------------------------------------------- */

//     std::vector <MKL_INT> perm (K_sc1.rows,0);
//     for (MKL_INT i = K_in.rows; i < K_sc1.rows; i++)
//     	perm[i] = 1;

//     MKL_INT nrhs = 0;

// //    phase = 11;
// //    PARDISO (pt, &maxfct, &mnum, &mtype, &phase,
// //        		&K_sc1.rows,
// //				&K_sc1.CSR_V_values[0], &K_sc1.CSR_I_row_indices[0], &K_sc1.CSR_J_col_indices[0],
// //				&perm[0], &nrhs, iparm, &msglvl, &ddum, &ddum, &error);
// //
// //    if ( error != 0 )
// //    {
// //    	printf ("\nERROR during symbolic factorization: %d", error);
// //    	exit(1);
// //    }


//     /* -------------------------------------------------------------------- */
//     /* .. Numerical factorization. */
//     /* -------------------------------------------------------------------- */

// 	SC_out.dense_values.resize(K_b_tmp.rows * K_b_tmp.rows);

//     phase = 12;
//     PARDISO (pt, &maxfct, &mnum, &mtype, &phase,
// 			&K_sc1.rows,
// 			&K_sc1.CSR_V_values[0], &K_sc1.CSR_I_row_indices[0], &K_sc1.CSR_J_col_indices[0],
// 			&perm[0], &nrhs,
// 			iparm, &msglvl, &ddum, &SC_out.dense_values[0], &error);

//     for (MKL_INT i = 0; i < SC_out.dense_values.size(); i++)
//     	SC_out.dense_values[i] = (-1.0)*SC_out.dense_values[i];

//     if ( error != 0 )
// 	{
// 		printf ("\nERROR during numerical factorization: %d", error);
// 		exit (2);
// 	} else {
// 		initialized = true;
// 	}

// 	/* -------------------------------------------------------------------- */
// 	/* .. Termination and release of memory. */
// 	/* -------------------------------------------------------------------- */
// 	phase = -1;           /* Release internal memory. */
// 	PARDISO (pt, &maxfct, &mnum, &mtype, &phase,
// 			&K_sc1.rows, &ddum, &K_sc1.CSR_I_row_indices[0], &K_sc1.CSR_J_col_indices[0], &idum, &nrhs,
// 			 iparm, &msglvl, &ddum, &ddum, &error);

// 	initialized = false;

//     /* -------------------------------------------------------------------- */
//     /* ..  allocate memory for the Schur-complement and copy it there.      */
//     /* -------------------------------------------------------------------- */

//     SC_out.cols = K_b_tmp.rows;
//     SC_out.rows = B0_in.cols;
//     SC_out.type = 'G';
//     SC_out.dense_values.resize( SC_out.cols * SC_out.rows);

//     //SC_out.ConvertDenseToCSR(1);

}


#include <stdio.h>
#include "mkl_rci.h"
#include "mkl_blas.h"
#include "mkl_spblas.h"
#include "mkl_service.h"

void SparseSolverCUDA::SolveCG(SparseMatrix & A_in, SEQ_VECTOR <double> & rhs_sol) {

	MKL_INT size = A_in.rows;
	SEQ_VECTOR <double> sol (size, 0);

	SolveCG(A_in, rhs_sol, sol);

	rhs_sol = sol;
}

void SparseSolverCUDA::SolveCG(SparseMatrix & A_in, SEQ_VECTOR <double> & rhs_in, SEQ_VECTOR <double> & sol) {
	SEQ_VECTOR<double> init;
	SolveCG(A_in, rhs_in, sol, init);
}

void SparseSolverCUDA::SolveCG(SparseMatrix & A_in, SEQ_VECTOR <double> & rhs_in, SEQ_VECTOR <double> & sol, SEQ_VECTOR <double> & initial_guess) {

	printf("Method SolveCG is not implemented yet.\n");
	exit(1);


// 	  /*---------------------------------------------------------------------------   */
// 	  /* Define arrays for the upper triangle of the coefficient matrix and rhs vector */
// 	  /* Compressed sparse row storage is used for sparse representation              */
// 	  /*---------------------------------------------------------------------------   */
// 	  MKL_INT rci_request, itercount, i;

// 	  //MKL_INT expected_itercount = 8;
// 	  MKL_INT n = A_in.rows; //rhs_in.size();


// 	  /* Fill all arrays containing matrix data. */
// 	  MKL_INT * ia  = &A_in.CSR_I_row_indices[0];
// 	  MKL_INT * ja  = &A_in.CSR_J_col_indices[0];
// 	  double  * a   = &A_in.CSR_V_values[0];


// //	  MKL_INT * ia  = CSR_I_row_indices;
// //	  MKL_INT * ja  = CSR_J_col_indices;
// //	  double  * a   = CSR_V_values;

// 	  /*---------------------------------------------------------------------------*/
// 	  /* Allocate storage for the solver ?par and temporary storage tmp            */
// 	  /*---------------------------------------------------------------------------*/
// 	  MKL_INT length = 128;
// 	  MKL_INT ipar[128];
// 	  double dpar[128];
// 	  char matdes[3];
// 	  double one = 1.E0;

// 	  SEQ_VECTOR <double> tmp_vec (4 * n, 0);
// 	  double * tmp;
// 	  tmp = &tmp_vec[0];

// 	  /*---------------------------------------------------------------------------*/
// 	  /* Some additional variables to use with the RCI (P)CG solver                */
// 	  /*---------------------------------------------------------------------------*/
// 	  double * solution = &sol[0];
//       double * rhs = &rhs_in[0];


// 	  double * temp;
// 	  SEQ_VECTOR <double> temp_vec (n,0);
// 	  temp = &temp_vec[0];

// 	  double euclidean_norm;

// 	  char tr = 'u';
// 	  double eone = -1.E0;
// 	  MKL_INT ione = 1;

// 	  /*---------------------------------------------------------------------------*/
// 	  /* Initialize the initial guess                                              */
// 	  /*---------------------------------------------------------------------------*/
// 	  if (initial_guess.size() > 0 ) {
// 		  for (i = 0; i < n; i++)
// 			solution[i] = initial_guess[i];
// 	  } else {
// 		  for (i = 0; i < n; i++)
// 			  solution[i] = 0.E0;
// 	  }
// 	  matdes[0] = 'd';
// 	  matdes[1] = 'l';
// 	  matdes[2] = 'n';
// 	  /*---------------------------------------------------------------------------*/
// 	  /* Initialize the solver                                                     */
// 	  /*---------------------------------------------------------------------------*/

// 	  dcg_init (&n, solution, rhs, &rci_request, ipar, dpar, tmp);
// 	  if (rci_request != 0)
// 	    goto failure;
// 	  /*---------------------------------------------------------------------------*/
// 	  /* Set the desired parameters:                                               */
// 	  /* LOGICAL parameters:                                                       */
// 	  /* -                                                                         */
// 	  /* INTEGER parameters:                                                       */
// 	  /* set the maximal number of iterations to 100                               */
// 	  /* DOUBLE parameters                                                         */
// 	  /* -                                                                         */
// 	  /*---------------------------------------------------------------------------*/
// 	  ipar[4] = 10000;
// 	  ipar[10]=1;
// 	  /*---------------------------------------------------------------------------*/
// 	  /* Check the correctness and consistency of the newly set parameters         */
// 	  /*---------------------------------------------------------------------------*/
// 	  dcg_check (&n, solution, rhs, &rci_request, ipar, dpar, tmp);
// 	  if (rci_request != 0)
// 	    goto failure;
// 	  /*---------------------------------------------------------------------------*/
// 	  /* Compute the solution by RCI (P)CG solver                                  */
// 	  /* Reverse Communications starts here                                        */
// 	  /*---------------------------------------------------------------------------*/
// 	rci:dcg (&n, solution, rhs, &rci_request, ipar, dpar, tmp);
// 	  ---------------------------------------------------------------------------
// 	  /* If rci_request=0, then the solution was found according to the requested  */
// 	  /* stopping tests. In this case, this means that it was found after 100      */
// 	  /* iterations.                                                               */
// 	  /*---------------------------------------------------------------------------*/
// 	  if (rci_request == 0)
// 	    goto getsln;
// 	  /*---------------------------------------------------------------------------*/
// 	  /* If rci_request=1, then compute the vector A*TMP[0]                        */
// 	  /* and put the result in vector TMP[n]                                       */
// 	  /*---------------------------------------------------------------------------*/
// 	  if (rci_request == 1)
// 	    {
// 	      mkl_dcsrsymv (&tr, &n, a, ia, ja, tmp, &tmp[n]);
// 	      goto rci;
// 	    }
// 	  /*---------------------------------------------------------------------------*/
// 	  /* If rci_request=2, then do the user-defined stopping test: compute the     */
// 	  /* Euclidean norm of the actual residual using MKL routines and check if     */
// 	  /* it is less than 1.E-8                                                     */
// 	  /*---------------------------------------------------------------------------*/
// 	  if (rci_request == 2)
// 	    {
// 	      mkl_dcsrsymv (&tr, &n, a, ia, ja, solution, temp);
// 	      daxpy (&n, &eone, rhs, &ione, temp, &ione);
// 	      euclidean_norm = dnrm2 (&n, temp, &ione);
// 	      /*---------------------------------------------------------------------------*/
// 	      /* The solution has not been found yet according to the user-defined stopping */
// 	      /* test. Continue RCI (P)CG iterations.                                      */
// 	      /*---------------------------------------------------------------------------*/
// 	      if (euclidean_norm > 1.E-13)
// 	        goto rci;
// 	      /*---------------------------------------------------------------------------*/
// 	      /* The solution has been found according to the user-defined stopping test   */
// 	      /*---------------------------------------------------------------------------*/
// 	      else
// 	        goto getsln;
// 	    }
// 	  /*---------------------------------------------------------------------------*/
// 	  /* If rci_request=3, then compute apply the preconditioner matrix C_inverse  */
// 	  /* on vector tmp[2*n] and put the result in vector tmp[3*n]                  */
// 	  /*---------------------------------------------------------------------------*/
// 	  if (rci_request == 3)
// 	    {
// 	      mkl_dcsrsv (&matdes[2], &n, &one, matdes, a, ja, ia, &ia[1], &tmp[2 * n], &tmp[3 * n]);
// 	      goto rci;
// 	    }	  /*---------------------------------------------------------------------------*/
// 	  /* If rci_request=anything else, then dcg subroutine failed                  */
// 	  /* to compute the solution vector: solution[n]                               */
// 	  /*---------------------------------------------------------------------------*/
// 	  goto failure;
// 	  /*---------------------------------------------------------------------------*/
// 	  /* Reverse Communication ends here                                           */
// 	  /* Get the current iteration number into itercount                           */
// 	  /*---------------------------------------------------------------------------*/
// 	getsln:dcg_get (&n, solution, rhs, &rci_request, ipar, dpar, tmp, &itercount);

// 	std::cout << " " << itercount;

// 	  //printf ("\nNumber of iterations: %d\n", itercount);

// 	  /*-------------------------------------------------------------------------*/
// 	  /* Release internal MKL memory that might be used for computations         */
// 	  /* NOTE: It is important to call the routine below to avoid memory leaks   */
// 	  /* unless you disable MKL Memory Manager                                   */
// 	  /*-------------------------------------------------------------------------*/
// 	failure: ; //printf ("This example FAILED as the solver has returned the ERROR code %d", rci_request);
// 	  //MKL_Free_Buffers ();


}
