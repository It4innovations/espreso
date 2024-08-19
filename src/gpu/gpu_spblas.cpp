
#if !defined(HAVE_CUDA) && !defined(HAVE_ROCM) && !defined(HAVE_ONEAPI)

#include "gpu_spblas.h"

namespace espreso {
namespace gpu {
namespace spblas {

    spblas_wrapper_impl get_implementation()
    {
        return spblas_wrapper_impl::NONE;
    }
    
    struct _handle {};

    struct _descr_matrix_csr {};

    struct _descr_matrix_dense {};

    struct _descr_vector_dense {};

    struct _descr_sparse_trsv {};

    struct _descr_sparse_trsm {};

    struct _descr_sparse_mv {};

    void handle_create(handle & h, mgm::queue & q) {}

    void handle_destroy(handle & h) {}

    template<typename T, typename I>
    void descr_matrix_csr_create(descr_matrix_csr & descr, I nrows, I ncols, I nnz, char fill) {}

    template<typename T, typename I, typename A>
    void descr_matrix_csr_link_data(descr_matrix_csr & descr, Matrix_CSR<T,I,A> & matrix) {}

    void descr_matrix_csr_destroy(descr_matrix_csr & descr) {}

    template<typename T, typename I>
    void descr_matrix_dense_create(descr_matrix_dense & descr, I nrows, I ncols, I ld, char order) {}

    template<typename T, typename I, typename A>
    void descr_matrix_dense_link_data(descr_matrix_dense & descr, Matrix_Dense<T,I,A> & matrix) {}

    void descr_matrix_dense_destroy(descr_matrix_dense & descr) {}

    template<typename T, typename I>
    void descr_vector_dense_create(descr_vector_dense & descr, I size) {}

    template<typename T, typename I, typename A>
    void descr_vector_dense_link_data(descr_vector_dense & descr, Vector_Dense<T,I,A> & vector) {}

    template<typename T, typename I, typename A>
    void descr_vector_dense_link_data(descr_vector_dense & descr, Matrix_Dense<T,I,A> & matrix, I colidx) {}

    void descr_vector_dense_destroy(descr_vector_dense & descr) {}

    void descr_sparse_trsv_create(descr_sparse_trsv & descr) {}

    void descr_sparse_trsv_destroy(descr_sparse_trsv & descr) {}

    void descr_sparse_trsm_create(descr_sparse_trsm & descr) {}

    void descr_sparse_trsm_destroy(descr_sparse_trsm & descr) {}

    void descr_sparse_mv_create(descr_sparse_mv & descr) {}

    void descr_sparse_mv_destroy(descr_sparse_mv & descr) {}

    template<typename T, typename I>
    void transpose(handle & h, descr_matrix_csr & output, descr_matrix_csr & input, bool conjugate, size_t & buffersize, void * buffer, char stage) {}

    template<typename T, typename I>
    void sparse_to_dense(handle & h, char transpose, descr_matrix_csr & sparse, descr_matrix_dense & dense, size_t & buffersize, void * buffer, char stage) {}

    template<typename T, typename I>
    void trsv(handle & h, char transpose, descr_matrix_csr & matrix, descr_vector_dense & rhs, descr_vector_dense & sol, descr_sparse_trsv & descr_trsv, size_t & buffersize, void * buffer, char stage) {}

    template<typename T, typename I>
    void trsm(handle & h, char transpose_mat, char transpose_rhs, char transpose_sol, descr_matrix_csr & matrix, descr_matrix_dense & rhs, descr_matrix_dense & sol, descr_sparse_trsm & descr_trsm, size_t & buffersize, void * buffer, char stage) {}

    template<typename T, typename I>
    void mv(handle & h, char transpose, descr_matrix_csr & A, descr_vector_dense & x, descr_vector_dense & y, descr_sparse_mv & descr_mv, size_t & buffersize, void * buffer, char stage) {}

    template<typename T, typename I>
    void mm(handle & h, char transpose_A, char transpose_B, descr_matrix_csr & A, descr_matrix_dense & B, descr_matrix_dense & C, size_t & buffersize, void * buffer, char stage) {}

}
}
}

#include "gpu/gpu_spblas.inst.hpp"

#endif
