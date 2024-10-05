
#ifndef SRC_GPU_SPBLAS_H_
#define SRC_GPU_SPBLAS_H_

#include "gpu_management.h"
#include "basis/utilities/cbmb_allocator.h"
#include "math/primitives/vector_dense.h"
#include "math/primitives/matrix_dense.h"
#include "math/primitives/matrix_csr.h"

namespace espreso {
namespace gpu {
namespace spblas {

    enum struct spblas_wrapper_impl {
        NONE,
        CUSPARSE_LEGACY,
        CUSPARSE_MODERN,
        ROCSPARSE,
        ONEMKL_SPARSE
    };

    // if function is in-place, both the in and out parameters have to be passed and be the same
    // if function is out-of-place, in and out parameters must not be the same
    enum struct place {
        NONE,
        IN_PLACE,
        OUT_OF_PLACE
    };

    spblas_wrapper_impl get_implementation();

    struct _handle;
    using handle = std::shared_ptr<_handle>;

    struct _descr_matrix_csr;
    using descr_matrix_csr = std::shared_ptr<_descr_matrix_csr>;

    struct _descr_matrix_dense;
    using descr_matrix_dense = std::shared_ptr<_descr_matrix_dense>;

    struct _descr_vector_dense;
    using descr_vector_dense = std::shared_ptr<_descr_vector_dense>;

    struct _descr_sparse_trsv;
    using descr_sparse_trsv = std::shared_ptr<_descr_sparse_trsv>;

    struct _descr_sparse_trsm;
    using descr_sparse_trsm = std::shared_ptr<_descr_sparse_trsm>;

    struct _descr_sparse_mv;
    using descr_sparse_mv = std::shared_ptr<_descr_sparse_mv>;

    struct buffer_sizes
    {
        size_t persistent = 0;
        size_t tmp_preprocess = 0;
        size_t tmp_update = 0;
        size_t tmp_compute = 0;
    };

    void handle_create(handle & h, mgm::queue & q);

    void handle_destroy(handle & h);

    template<typename T, typename I>
    void descr_matrix_csr_create(handle & h, descr_matrix_csr & descr, I nrows, I ncols, I nnz, char fill);

    template<typename T, typename I, typename A>
    void descr_matrix_csr_link_data(handle & h, descr_matrix_csr & descr, Matrix_CSR<T,I,A> & matrix);

    void descr_matrix_csr_destroy(handle & h, descr_matrix_csr & descr);

    template<typename T, typename I>
    void descr_matrix_dense_create(handle & h, descr_matrix_dense & descr, I nrows, I ncols, I ld, char order);

    template<typename T, typename I, typename A>
    void descr_matrix_dense_link_data(handle & h, descr_matrix_dense & descr, Matrix_Dense<T,I,A> & matrix);

    void descr_matrix_dense_destroy(handle & h, descr_matrix_dense & descr);

    template<typename T, typename I>
    void descr_vector_dense_create(handle & h, descr_vector_dense & descr, I size);

    template<typename T, typename I, typename A>
    void descr_vector_dense_link_data(handle & h, descr_vector_dense & descr, Vector_Dense<T,I,A> & vector);

    template<typename T, typename I, typename A>
    void descr_vector_dense_link_data(handle & h, descr_vector_dense & descr, Matrix_Dense<T,I,A> & matrix, I colidx = 0);

    void descr_vector_dense_destroy(handle & h, descr_vector_dense & descr);

    void descr_sparse_trsv_create(handle & h, descr_sparse_trsv & descr);

    void descr_sparse_trsv_destroy(handle & h, descr_sparse_trsv & descr);

    void descr_sparse_trsm_create(handle & h, descr_sparse_trsm & descr);

    void descr_sparse_trsm_destroy(handle & h, descr_sparse_trsm & descr);

    void descr_sparse_mv_create(handle & h, descr_sparse_mv & descr);

    void descr_sparse_mv_destroy(handle & h, descr_sparse_mv & descr);

    place get_place_trsm();

    // stages: Buffersize, Preprocess, Compute
    template<typename T, typename I>
    void transpose(handle & h, descr_matrix_csr & output, descr_matrix_csr & input, bool conjugate, size_t & buffersize, void * buffer, char stage);

    // stages: Buffersize, Compute
    template<typename T, typename I>
    void sparse_to_dense(handle & h, char transpose, descr_matrix_csr & sparse, descr_matrix_dense & dense, size_t & buffersize, void * buffer, char stage);

    // stages: Buffersize, Preprocess, Update, Compute
    template<typename T, typename I>
    void trsv(handle & h, char transpose, descr_matrix_csr & matrix, descr_vector_dense & rhs, descr_vector_dense & sol, descr_sparse_trsv & descr_trsv, buffer_sizes & buffersizes, void * buffer_persistent, void * buffer_tmp, char stage);

    // stages: Buffersize, Preprocess, Update, Compute
    template<typename T, typename I>
    void trsm(handle & h, char transpose_mat, char transpose_rhs, char transpose_sol, descr_matrix_csr & matrix, descr_matrix_dense & rhs, descr_matrix_dense & sol, descr_sparse_trsm & descr_trsm, buffer_sizes & buffersizes, void * buffer_persistent, void * buffer_tmp, char stage);

    // stages: Buffersize, Preprocess, Compute
    template<typename T, typename I>
    void mv(handle & h, char transpose, descr_matrix_csr & A, descr_vector_dense & x, descr_vector_dense & y, descr_sparse_mv & descr_mv, size_t & buffersize, void * buffer, char stage);

    // stages: Buffersize, Preprocess, Compute
    template<typename T, typename I>
    void mm(handle & h, char transpose_A, char transpose_B, descr_matrix_csr & A, descr_matrix_dense & B, descr_matrix_dense & C, size_t & buffersize, void * buffer, char stage);

}
}
}

#endif /* SRC_GPU_SPBLAS_H_ */
