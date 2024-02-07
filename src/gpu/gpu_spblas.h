
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

    struct _handle;
    struct handle { _handle *inner; };

    struct _descr_matrix_csr;
    struct _descr_matrix_dense;
    struct _descr_vector_dense;
    struct _descr_sparse_trsv;
    struct _descr_sparse_trsm;

    struct descr_matrix_csr   { _descr_matrix_csr   *inner; };
    struct descr_matrix_dense { _descr_matrix_dense *inner; };
    struct descr_vector_dense { _descr_vector_dense *inner; };
    struct descr_sparse_trsv  { _descr_sparse_trsv  *inner; };
    struct descr_sparse_trsm  { _descr_sparse_trsm  *inner; };

    void handle_create(handle & h, mgm::queue & q);

    void handle_destroy(handle & h);

    template<typename T, typename I>
    void descr_matrix_csr_create(descr_matrix_csr & descr, I nrows, I ncols, I nnz, char symmetry);

    template<typename T, typename I, typename A>
    void descr_matrix_sparse_link_data(descr_matrix_csr & descr, Matrix_CSR<T,I,A> & matrix);

    void descr_matrix_csr_destroy(descr_matrix_csr & descr);

    template<typename T, typename I>
    void descr_matrix_dense_create(descr_matrix_dense & descr, I nrows, I ncols, I ld, char order);

    template<typename T, typename I, typename A>
    void descr_matrix_dense_link_data(descr_matrix_dense & descr, Matrix_Dense<T,I,A> & matrix);

    void descr_matrix_dense_destroy(descr_matrix_dense & descr);

    template<typename T, typename I>
    void descr_vector_dense_create(descr_vector_dense & descr, I size);

    template<typename T, typename I, typename A>
    void descr_vector_dense_link_data(descr_vector_dense & descr, Vector_Dense<T,I,A> & vector);

    template<typename T, typename I, typename A>
    void descr_vector_dense_link_data(descr_vector_dense & descr, Matrix_Dense<T,I,A> & matrix, I colidx = 0);

    void descr_vector_dense_destroy(descr_vector_dense & descr);

    void descr_sparse_trsv_create(descr_sparse_trsv & descr);

    void descr_sparse_trsv_destroy(descr_sparse_trsv & descr);

    void descr_sparse_trsm_create(descr_sparse_trsm & descr);

    void descr_sparse_trsm_destroy(descr_sparse_trsm & descr);

    template<typename T, typename I>
    void sparse_to_dense(handle & h, char transpose, descr_matrix_csr & sparse, descr_matrix_dense & dense, size_t & buffersize, void * buffer, char stage);

    template<typename T, typename I>
    void trsv(handle & h, char transpose, descr_matrix_csr & matrix, descr_vector_dense & rhs, descr_vector_dense & sol, descr_sparse_trsv & descr_trsv, size_t & buffersize, void * buffer, char stage);

    template<typename T, typename I>
    void trsm(handle & h, char transpose_mat, char transpose_rhs, descr_matrix_csr & matrix, descr_matrix_dense & rhs, descr_matrix_dense & sol, descr_sparse_trsm & descr_trsm, size_t & buffersize, void * buffer, char stage);

    template<typename T, typename I>
    static void mv(handle & h, char transpose, descr_matrix_csr & A, descr_vector_dense & x, descr_vector_dense & y, size_t & buffersize, void * buffer, char stage);

    template<typename T, typename I>
    void mm(handle & h, char transpose_A, char transpose_B, descr_matrix_csr & A, descr_matrix_dense & B, descr_matrix_dense & C, size_t & buffersize, void * buffer, char stage);

}
}
}

#endif /* SRC_GPU_SPBLAS_H_ */
