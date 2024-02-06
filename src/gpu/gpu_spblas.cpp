
#include "gpu_spblas.h"

namespace espreso {
namespace gpu {
namespace spblas {

#ifndef HAVE_CUDA

void handle_create(handle & h, mgm::queue & q) {}

void handle_destroy(handle & h) {}

template<>
void descr_matrix_csr_create<double, int>(descr_matrix_csr & descr, int nrows, int ncols, int nnz, char symmetry) {}

template<>
void descr_matrix_sparse_link_data<double, int, mgm::Ah>(descr_matrix_csr & descr, Matrix_CSR<double,int,mgm::Ah> & matrix) {}

template<>
void descr_matrix_sparse_link_data<double, int, mgm::Ad>(descr_matrix_csr & descr, Matrix_CSR<double,int,mgm::Ad> & matrix) {}

void descr_matrix_csr_destroy(descr_matrix_csr & descr) {}

template<>
void descr_matrix_dense_create<double, int>(descr_matrix_dense & descr, int nrows, int ncols, int ld, char order) {}

template<>
void descr_matrix_dense_link_data(descr_matrix_dense & descr, Matrix_Dense<double,int,mgm::Ah> & matrix) {}

template<>
void descr_matrix_dense_link_data(descr_matrix_dense & descr, Matrix_Dense<double,int,mgm::Ad> & matrix) {}

template<>
void descr_matrix_dense_link_data(descr_matrix_dense & descr, Matrix_Dense<double,int,cbmba<true, false> > & matrix) {}

void descr_matrix_dense_destroy(descr_matrix_dense & descr) {}

template<>
void descr_vector_dense_create<double, int>(descr_vector_dense & descr, int size) {}

template<>
void descr_vector_dense_link_data(descr_vector_dense & descr, Vector_Dense<double,int,mgm::Ah> & vector) {}

template<>
void descr_vector_dense_link_data(descr_vector_dense & descr, Vector_Dense<double,int,mgm::Ad> & vector) {}

template<>
void descr_vector_dense_link_data(descr_vector_dense & descr, Vector_Dense<double,int,cbmba<true, false> > & vector) {}

template<>
void descr_vector_dense_link_data(descr_vector_dense & descr, Matrix_Dense<double,int,mgm::Ah> & matrix, int colidx) {}

template<>
void descr_vector_dense_link_data(descr_vector_dense & descr, Matrix_Dense<double,int,mgm::Ad> & matrix, int colidx) {}

template<>
void descr_vector_dense_link_data(descr_vector_dense & descr, Matrix_Dense<double,int,cbmba<true, false> > & matrix, int colidx) {}

void descr_vector_dense_destroy(descr_vector_dense & descr) {}

void descr_sparse_trsv_create(descr_sparse_trsv & descr) {}

void descr_sparse_trsv_destroy(descr_sparse_trsv & descr) {}

void descr_sparse_trsm_create(descr_sparse_trsm & descr) {}

void descr_sparse_trsm_destroy(descr_sparse_trsm & descr) {}

template<>
void sparse_to_dense<double, int>(handle & h, char transpose, descr_matrix_csr & sparse, descr_matrix_dense & dense, size_t & buffersize, void * buffer, char stage) {}

template<>
void trsv<double, int>(handle & h, char transpose, descr_matrix_csr & matrix, descr_vector_dense & rhs, descr_vector_dense & sol, descr_sparse_trsv & descr_trsv, size_t & buffersize, void * buffer, char stage) {}

template<>
void trsm<double, int>(handle & h, char transpose_mat, char transpose_rhs, descr_matrix_csr & matrix, descr_matrix_dense & rhs, descr_matrix_dense & sol, descr_sparse_trsm & descr_trsm, size_t & buffersize, void * buffer, char stage) {}

template<>
void mv<double, int>(handle & h, char transpose, descr_matrix_csr & A, descr_vector_dense & x, descr_vector_dense & y, size_t & buffersize, void * buffer, char stage) {}

template<>
void mm<double, int>(handle & h, char transpose_A, char transpose_B, descr_matrix_csr & A, descr_matrix_dense & B, descr_matrix_dense & C, size_t & buffersize, void * buffer, char stage) {}

#endif

}
}
}
