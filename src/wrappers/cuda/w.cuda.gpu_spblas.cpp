
#ifdef HAVE_CUDA

#ifdef USE_CUSPARSE_LEGACY
#include "w.cuda.gpu_spblas_legacy.hpp"
#else
#include "w.cuda.gpu_spblas_modern.hpp"
#endif

namespace espreso {
namespace gpu {
namespace spblas {

void handle_create(handle & h, mgm::queue & q)
{
    h.inner = new _handle();
    handle_create(*h.inner, *q.inner);
}

void handle_destroy(handle & h)
{
    handle_destroy(*h.inner);
    delete h.inner;
}

template<>
void descr_matrix_csr_create<double, int>(descr_matrix_csr & descr, int nrows, int ncols, int nnz, char symmetry)
{
    descr.inner = new _descr_matrix_csr();
    descr_matrix_csr_create<double>(*descr.inner, nrows, ncols, nnz, symmetry);
}

template<>
void descr_matrix_sparse_link_data<double, int, mgm::Ah>(descr_matrix_csr & descr, Matrix_CSR<double,int,mgm::Ah> & matrix)
{
    descr_matrix_sparse_link_data(*descr.inner, matrix);
}

template<>
void descr_matrix_sparse_link_data<double, int, mgm::Ad>(descr_matrix_csr & descr, Matrix_CSR<double,int,mgm::Ad> & matrix)
{
    descr_matrix_sparse_link_data(*descr.inner, matrix);
}

void descr_matrix_csr_destroy(descr_matrix_csr & descr)
{
    descr_matrix_csr_destroy(*descr.inner);
    delete descr.inner;
}

template<>
void descr_matrix_dense_create<double, int>(descr_matrix_dense & descr, int nrows, int ncols, int ld, char order)
{
    descr.inner = new _descr_matrix_dense();
    descr_matrix_dense_create<double, int>(*descr.inner, nrows, ncols, ld, order);
}

template<>
void descr_matrix_dense_link_data(descr_matrix_dense & descr, Matrix_Dense<double,int,mgm::Ah> & matrix)
{
    descr_matrix_dense_link_data(*descr.inner, matrix);
}

template<>
void descr_matrix_dense_link_data(descr_matrix_dense & descr, Matrix_Dense<double,int,mgm::Ad> & matrix)
{
    descr_matrix_dense_link_data(*descr.inner, matrix);
}

template<>
void descr_matrix_dense_link_data(descr_matrix_dense & descr, Matrix_Dense<double,int,cbmba<true, false> > & matrix)
{
    descr_matrix_dense_link_data(*descr.inner, matrix);
}

void descr_matrix_dense_destroy(descr_matrix_dense & descr)
{
    descr_matrix_dense_destroy(*descr.inner);
    delete descr.inner;
}

template<>
void descr_vector_dense_create<double, int>(descr_vector_dense & descr, int size)
{
    descr.inner = new _descr_vector_dense();
    descr_vector_dense_create<double, int>(*descr.inner, size);
}

template<>
void descr_vector_dense_link_data(descr_vector_dense & descr, Vector_Dense<double,int,mgm::Ah> & vector)
{
    descr_vector_dense_link_data(*descr.inner, vector);
}

template<>
void descr_vector_dense_link_data(descr_vector_dense & descr, Vector_Dense<double,int,mgm::Ad> & vector)
{
    descr_vector_dense_link_data(*descr.inner, vector);
}

template<>
void descr_vector_dense_link_data(descr_vector_dense & descr, Vector_Dense<double,int,cbmba<true, false> > & vector)
{
    descr_vector_dense_link_data(*descr.inner, vector);
}

template<>
void descr_vector_dense_link_data(descr_vector_dense & descr, Matrix_Dense<double,int,mgm::Ah> & matrix, int colidx)
{
    descr_vector_dense_link_data(*descr.inner, matrix, colidx);
}

template<>
void descr_vector_dense_link_data(descr_vector_dense & descr, Matrix_Dense<double,int,mgm::Ad> & matrix, int colidx)
{
    descr_vector_dense_link_data(*descr.inner, matrix, colidx);
}

template<>
void descr_vector_dense_link_data(descr_vector_dense & descr, Matrix_Dense<double,int,cbmba<true, false> > & matrix, int colidx)
{
    descr_vector_dense_link_data(*descr.inner, matrix, colidx);
}

void descr_vector_dense_destroy(descr_vector_dense & descr)
{
    descr_vector_dense_destroy(*descr.inner);
    delete descr.inner;
}

void descr_sparse_trsv_create(descr_sparse_trsv & descr)
{
    descr.inner = new _descr_sparse_trsv();
    descr_sparse_trsv_create(*descr.inner);
}

void descr_sparse_trsv_destroy(descr_sparse_trsv & descr)
{
    descr_sparse_trsv_destroy(*descr.inner);
    delete descr.inner;
}

void descr_sparse_trsm_create(descr_sparse_trsm & descr)
{
    descr.inner = new _descr_sparse_trsm();
    descr_sparse_trsm_create(*descr.inner);
}

void descr_sparse_trsm_destroy(descr_sparse_trsm & descr)
{
    descr_sparse_trsm_destroy(*descr.inner);
    delete descr.inner;
}

template<>
void sparse_to_dense<double, int>(handle & h, char transpose, descr_matrix_csr & sparse, descr_matrix_dense & dense, size_t & buffersize, void * buffer, char stage)
{
    sparse_to_dense<double, int>(*h.inner, transpose, *sparse.inner, *dense.inner, buffersize, buffer, stage);
}

template<>
void trsv<double, int>(handle & h, char transpose, descr_matrix_csr & matrix, descr_vector_dense & rhs, descr_vector_dense & sol, descr_sparse_trsv & descr_trsv, size_t & buffersize, void * buffer, char stage)
{
    trsv<double, int>(*h.inner, transpose, *matrix.inner, *rhs.inner, *sol.inner, *descr_trsv.inner, buffersize, buffer, stage);
}

template<>
void trsm<double, int>(handle & h, char transpose_mat, char transpose_rhs, descr_matrix_csr & matrix, descr_matrix_dense & rhs, descr_matrix_dense & sol, descr_sparse_trsm & descr_trsm, size_t & buffersize, void * buffer, char stage)
{
    trsm<double, int>(*h.inner, transpose_mat, transpose_rhs, *matrix.inner, *rhs.inner, *sol.inner, *descr_trsm.inner, buffersize, buffer, stage);
}

template<>
void mv<double, int>(handle & h, char transpose, descr_matrix_csr & A, descr_vector_dense & x, descr_vector_dense & y, size_t & buffersize, void * buffer, char stage)
{
    mv<double, int>(*h.inner, transpose, *A.inner, *x.inner, *y.inner, buffersize, buffer, stage);
}

template<>
void mm<double, int>(handle & h, char transpose_A, char transpose_B, descr_matrix_csr & A, descr_matrix_dense & B, descr_matrix_dense & C, size_t & buffersize, void * buffer, char stage)
{
    mm<double, int>(*h.inner, transpose_A, transpose_B, *A.inner, *B.inner, *C.inner, buffersize, buffer, stage);
}

}
}
}

#endif

