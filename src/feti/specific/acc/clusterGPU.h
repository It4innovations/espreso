/*
 * clusterGPU.h
 *
 *  Created on: Feb 24, 2016
 *      Author: lriha
 */

#ifndef SOLVER_SPECIFIC_ACC_CLUSTERGPU_H_
#define SOLVER_SPECIFIC_ACC_CLUSTERGPU_H_

#include "feti/specific/cluster.h"
// TODO: Should be moved to the CUDA wrapper:
#include <cuda.h>
#include <cuda_runtime.h>
#include "cusparse_v2.h"

namespace espreso {
class ClusterGPU: public ClusterBase
{

public:
	// Constructor
	ClusterGPU(const FETIConfiguration &configuration, DataHolder *instance_in): ClusterBase(configuration, instance_in), device_id(-1) { }
	~ClusterGPU();

	void Create_SC_perDomain(bool USE_FLOAT);
    void CreateDirichletPrec(DataHolder *instance);
	void SetupKsolvers();

	void multKplusGlobal_GPU   ( SEQ_VECTOR<SEQ_VECTOR<double> > & x_in );

	void multKplus_HF      (SEQ_VECTOR<SEQ_VECTOR<double> > & x_in);
	void multKplus_HF_SC   (SEQ_VECTOR<SEQ_VECTOR<double> > & x_in);
	void multKplus_HF_SC   (SEQ_VECTOR<SEQ_VECTOR<double> > & x_in, SEQ_VECTOR<SEQ_VECTOR<double> > & y_out);
	void multKplus_HF_SPDS (SEQ_VECTOR<SEQ_VECTOR<double> > & x_in);

	void multKplus_HF_Loop1 (SEQ_VECTOR<SEQ_VECTOR<double> > & x_in);
	void multKplus_HF_CP    ();

	void multKplus_HF_Loop2_SC   (SEQ_VECTOR<SEQ_VECTOR<double> > & x_in, SEQ_VECTOR<SEQ_VECTOR<double> > & y_out);
	void multKplus_HF_Loop2_SPDS (SEQ_VECTOR<SEQ_VECTOR<double> > & x_in);
	void multKplus_HF_Loop2_MIX  (SEQ_VECTOR<SEQ_VECTOR<double> > & x_in);

	void CreateCudaStreamPool();
	void DestroyCudaStreamPool();

	SEQ_VECTOR <cudaStream_t> cuda_stream_pool;
    SEQ_VECTOR <cublasHandle_t> cublas_handle_pool;

private:
	void GetSchurComplement(bool USE_FLOAT, esint i);
	// The method is prepared for multiple GPUs per cluster
	void GetSchurComplementsGpu(bool USE_FLOAT, SEQ_VECTOR<int>& vec_L_nnz,
     SEQ_VECTOR<int*>& vec_L_row_indexes, SEQ_VECTOR<int*>& vec_L_col_pointers, SEQ_VECTOR<double*>& vec_L_values,
     SEQ_VECTOR<int>& vec_U_nnz, SEQ_VECTOR<int*>& vec_U_row_indexes, SEQ_VECTOR<int*>& vec_U_col_pointers,
     SEQ_VECTOR<double*>& vec_U_values, SEQ_VECTOR<int*>& vec_perm, SEQ_VECTOR<int*>& vec_perm_2, esint max_B1_nnz,
     esint max_B1_rows, esint max_B1_size, esint max_K_rows, esint max_L_nnz, esint max_U_nnz);
    void GetSchurComplementsCpu(bool USE_FLOAT);
    void GetDirichletPrec(DataHolder *instance, esint d);
	void GetGPU();
    // Calculates GPU buffers apaort from dense LSC matrix = B1_comp_dom.rows * B1_comp_dom.rows * sizeof(double)
    // and in/out vectors = 2 * vec_size * sizeof(double)
    size_t CalculateGpuBufferSize(esint max_B1_nnz, esint max_B1_rows, esint max_B1_size, esint max_K_rows, esint max_L_nnz, esint max_U_nnz);
	// TODO change to arrays for multi-GPU per cluster
	size_t  GPU_free_mem;
	size_t  GPU_total_mem;
	int device_id;
	SEQ_VECTOR<int> lsc_on_gpu_ids;
	SEQ_VECTOR<int> lsc_on_cpu_ids;
    // numbeŕ of computation streams for LSC
    int n_streams_per_gpu = 2;
    // Determines how often the GPU is synchronized in order to free temporary memory
    int n_csrsm2_info_per_gpu; // Set in ecf
    // The struct represents members assigned to one GPU
    typedef struct {
        // Number of LSCs stored in each GPU
        int n_lsc_gpu;
        // Define ranges
        int start;
        int end;

        // 1 data_stream per GPU
        cudaStream_t data_stream;

        cudaStream_t *h_array_stream;
        cusparseHandle_t *h_array_handle;

        // Should be shared with itersolverGPU in the new dev branch
        cublasHandle_t cublas_handle;

        cudaEvent_t event_data_preload;
        cudaEvent_t event1;
        cudaEvent_t event2;

        csrsm2Info_t *h_array_info_L;
        csrsm2Info_t *h_array_info_U;

        // Only one set of device arrays for CSC L and U (L^T)
        int *d_csc_L_col_ptr;
        int *d_csc_L_row_ind;
        double *d_csc_L_val;

        int *d_csc_U_col_ptr;
        int *d_csc_U_row_ind;
        double *d_csc_U_val;

        double *d_sym_full_lsc;

        int **h_array_d_csr_B_row_ptr;
        int **h_array_d_csr_B_col_ind;
        double **h_array_d_csr_B_val;
        char **h_array_d_buffer;

        // Possible future micro-optimization: 
        // Comment out in order to eliminate Bt_csr (this would currently break HtoD data transfers)
        int **h_array_d_csr_Bt_row_ptr;
        int **h_array_d_csr_Bt_col_ind;
        double **h_array_d_csr_Bt_val;

        double **h_array_d_X_reordered;
        double **h_array_d_lsc;
        int *h_array_lsc_id;

        // Permutation vectors
        int **h_array_d_pinv;
        int **h_array_d_q;
    } TGPU;

// TODO Distribution of domains based on real size of LSC of each domain (# RHS)
// e.g.: Having N domains and G GPUs, sort domains descending by LSC size
// assign G largest domains to GPUs, and continue assigning the following domains
// to the GPU with the currently smallest sum of assigned LSC sizes.
//  add int* lsc_sizes argument
void DistributeDomains(TGPU* gpus, int n_gpu, int n_lsc);

#ifdef SHARE_SC
	SEQ_VECTOR <double *> SC_dense_val_orig;
	SEQ_VECTOR <float *> SC_dense_val_orig_fl;
#endif
};

}

#endif /* SOLVER_SPECIFIC_ACC_CLUSTERGPU_H_ */
