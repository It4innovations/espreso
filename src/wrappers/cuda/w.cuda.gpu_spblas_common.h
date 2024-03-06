
#ifndef SRC_WRAPPERS_CUDA_W_CUDA_GPU_SPBLAS_COMMON_H
#define SRC_WRAPPERS_CUDA_W_CUDA_GPU_SPBLAS_COMMON_H



namespace espreso {
namespace gpu {
namespace spblas {

    template<typename I>
    void my_csr_transpose_buffersize(cudaStream_t & stream, I input_nrows, I input_ncols, I nnz, size_t & buffersize);
    
    template<typename I>
    void my_csr_transpose_preprocess(cudaStream_t & stream, I input_nrows, I input_ncols, I nnz, I * input_rowptrs, I * input_colidxs, I * output_rowptrs, I * output_colidxs, size_t buffersize, void * buffer);
    
    template<typename T, typename I>
    void my_csr_transpose_compute(cudaStream_t & stream, I nnz, T * input_vals, T * output_vals, void * buffer);

}
}
}

#endif
