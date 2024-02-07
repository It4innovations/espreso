
#if !defined(HAVE_CUDA) && !defined(HAVE_ROCM)

#include "gpu_spblas.h"

namespace espreso {
namespace gpu {
namespace spblas {
    
    struct _handle {};

    struct _descr_matrix_csr {};

    struct _descr_matrix_dense {};

    struct _descr_vector_dense {};

    struct _descr_sparse_trsv {};

    struct _descr_sparse_trsm {};

    void handle_create(handle & h, mgm::queue & q) {}

    void handle_destroy(handle & h) {}

    template<typename T, typename I>
    void descr_matrix_csr_create(descr_matrix_csr & descr, I nrows, I ncols, I nnz, char symmetry) {}

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
    void descr_vector_dense_link_data(descr_vector_dense & descr, Matrix_Dense<T,I,A> & matrix, I colidx = 0) {}

    void descr_vector_dense_destroy(descr_vector_dense & descr) {}

    void descr_sparse_trsv_create(descr_sparse_trsv & descr) {}

    void descr_sparse_trsv_destroy(descr_sparse_trsv & descr) {}

    void descr_sparse_trsm_create(descr_sparse_trsm & descr) {}

    void descr_sparse_trsm_destroy(descr_sparse_trsm & descr) {}

    template<typename T, typename I>
    void sparse_to_dense(handle & h, char transpose, descr_matrix_csr & sparse, descr_matrix_dense & dense, size_t & buffersize, void * buffer, char stage) {}

    template<typename T, typename I>
    void trsv(handle & h, char transpose, descr_matrix_csr & matrix, descr_vector_dense & rhs, descr_vector_dense & sol, descr_sparse_trsv & descr_trsv, size_t & buffersize, void * buffer, char stage) {}

    template<typename T, typename I>
    void trsm(handle & h, char transpose_mat, char transpose_rhs, descr_matrix_csr & matrix, descr_matrix_dense & rhs, descr_matrix_dense & sol, descr_sparse_trsm & descr_trsm, size_t & buffersize, void * buffer, char stage) {}

    template<typename T, typename I>
    void mv(handle & h, char transpose, descr_matrix_csr & A, descr_vector_dense & x, descr_vector_dense & y, size_t & buffersize, void * buffer, char stage) {}

    template<typename T, typename I>
    void mm(handle & h, char transpose_A, char transpose_B, descr_matrix_csr & A, descr_matrix_dense & B, descr_matrix_dense & C, size_t & buffersize, void * buffer, char stage) {}



    #define INSTANTIATE(T,I,Adevice) \
    template void descr_matrix_csr_create<T,I>(descr_matrix_csr & descr, I nrows, I ncols, I nnz, char symmetry); \
    template void descr_matrix_csr_link_data<T,I,Adevice>(descr_matrix_csr & descr, Matrix_CSR<T,I,A> & matrix); \
    template void descr_matrix_dense_create<T,I>(descr_matrix_dense & descr, I nrows, I ncols, I ld, char order); \
    template void descr_matrix_dense_link_data<T,I,Adevice>(descr_matrix_dense & descr, Matrix_Dense<T,I,A> & matrix); \
    template void descr_vector_dense_create<T,I>(descr_vector_dense & descr, I size); \
    template void descr_vector_dense_link_data<T,I,Adevice>(descr_vector_dense & descr, Vector_Dense<T,I,A> & vector); \
    template void descr_vector_dense_link_data<T,I,Adevice>(descr_vector_dense & descr, Matrix_Dense<T,I,A> & matrix, I colidx = 0); \
    template void sparse_to_dense<T,I>(handle & h, char transpose, descr_matrix_csr & sparse, descr_matrix_dense & dense, size_t & buffersize, void * buffer, char stage); \
    template void trsv<T,I>(handle & h, char transpose, descr_matrix_csr & matrix, descr_vector_dense & rhs, descr_vector_dense & sol, descr_sparse_trsv & descr_trsv, size_t & buffersize, void * buffer, char stage); \
    template void trsm<T,I>(handle & h, char transpose_mat, char transpose_rhs, descr_matrix_csr & matrix, descr_matrix_dense & rhs, descr_matrix_dense & sol, descr_sparse_trsm & descr_trsm, size_t & buffersize, void * buffer, char stage); \
    template void mv<T,I>(handle & h, char transpose, descr_matrix_csr & A, descr_vector_dense & x, descr_vector_dense & y, size_t & buffersize, void * buffer, char stage); \
    template void mm<T,I>(handle & h, char transpose_A, char transpose_B, descr_matrix_csr & A, descr_matrix_dense & B, descr_matrix_dense & C, size_t & buffersize, void * buffer, char stage);
        // INSTANTIATE(float,                int32_t, mgm::Ad)
        INSTANTIATE(double,               int32_t, mgm::Ad)
        // INSTANTIATE(std::complex<float >, int32_t, mgm::Ad)
        // INSTANTIATE(std::complex<double>, int32_t, mgm::Ad)
        // INSTANTIATE(float,                int64_t, mgm::Ad)
        // INSTANTIATE(double,               int64_t, mgm::Ad)
        // INSTANTIATE(std::complex<float >, int64_t, mgm::Ad)
        // INSTANTIATE(std::complex<double>, int64_t, mgm::Ad)
        // INSTANTIATE(float,                int32_t, cbmba_d)
        INSTANTIATE(double,               int32_t, cbmba_d)
        // INSTANTIATE(std::complex<float >, int32_t, cbmba_d)
        // INSTANTIATE(std::complex<double>, int32_t, cbmba_d)
        // INSTANTIATE(float,                int64_t, cbmba_d)
        // INSTANTIATE(double,               int64_t, cbmba_d)
        // INSTANTIATE(std::complex<float >, int64_t, cbmba_d)
        // INSTANTIATE(std::complex<double>, int64_t, cbmba_d)
    #undef INSTANTIATE


}
}
}

#endif
