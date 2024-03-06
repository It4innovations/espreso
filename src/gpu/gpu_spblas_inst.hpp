
namespace espreso {
namespace gpu {
namespace spblas {

    #define INSTANTIATE_T_I_ADEVICE(T,I,Adevice) \
    template void descr_matrix_csr_link_data<T,I,Adevice>(descr_matrix_csr & descr, Matrix_CSR<T,I,Adevice> & matrix); \
    template void descr_matrix_dense_link_data<T,I,Adevice>(descr_matrix_dense & descr, Matrix_Dense<T,I,Adevice> & matrix); \
    template void descr_vector_dense_link_data<T,I,Adevice>(descr_vector_dense & descr, Vector_Dense<T,I,Adevice> & vector); \
    template void descr_vector_dense_link_data<T,I,Adevice>(descr_vector_dense & descr, Matrix_Dense<T,I,Adevice> & matrix, I colidx = 0);

        #define INSTANTIATE_T_I(T,I) \
        template void descr_matrix_csr_create<T,I>(descr_matrix_csr & descr, I nrows, I ncols, I nnz, char fill); \
        template void descr_matrix_dense_create<T,I>(descr_matrix_dense & descr, I nrows, I ncols, I ld, char order); \
        template void descr_vector_dense_create<T,I>(descr_vector_dense & descr, I size); \
        template void transpose<T,I>(handle & h, descr_matrix_csr & output, descr_matrix_csr & input, size_t & buffersize, void * buffer, char stage); \
        template void sparse_to_dense<T,I>(handle & h, char transpose, descr_matrix_csr & sparse, descr_matrix_dense & dense, size_t & buffersize, void * buffer, char stage); \
        template void trsv<T,I>(handle & h, char transpose, descr_matrix_csr & matrix, descr_vector_dense & rhs, descr_vector_dense & sol, descr_sparse_trsv & descr_trsv, size_t & buffersize, void * buffer, char stage); \
        template void trsm<T,I>(handle & h, char transpose_mat, char transpose_rhs, char transpose_sol, descr_matrix_csr & matrix, descr_matrix_dense & rhs, descr_matrix_dense & sol, descr_sparse_trsm & descr_trsm, size_t & buffersize, void * buffer, char stage); \
        template void mv<T,I>(handle & h, char transpose, descr_matrix_csr & A, descr_vector_dense & x, descr_vector_dense & y, size_t & buffersize, void * buffer, char stage); \
        template void mm<T,I>(handle & h, char transpose_A, char transpose_B, descr_matrix_csr & A, descr_matrix_dense & B, descr_matrix_dense & C, size_t & buffersize, void * buffer, char stage); \
        INSTANTIATE_T_I_ADEVICE(T, I, mgm::Ad) \
        INSTANTIATE_T_I_ADEVICE(T, I, cbmba_d)

            #define INSTANTIATE_T(T) \
            INSTANTIATE_T_I(T, int32_t) \
            /* INSTANTIATE_T_I(T, int64_t) */

                // INSTANTIATE_T(float)
                INSTANTIATE_T(double)
                // INSTANTIATE_T(std::complex<float>)
                // INSTANTIATE_T(std::complex<double>)

            #undef INSTANTIATE_T
        #undef INSTANTIATE_T_I
    #undef INSTANTIATE_T_I_ADEVICE

}
}
}
