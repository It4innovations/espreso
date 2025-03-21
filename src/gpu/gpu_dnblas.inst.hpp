
namespace espreso {
namespace gpu {
namespace dnblas {

    #define INSTANTIATE_T_I(T,I) \
    template void trsv<T,I>(handle & h, I n, T * A, I ld_A, char order_A, char op_A, char fill_A, T * x); \
    template void trsm<T,I>(handle & h, char side, I n, I nrhs, T * A, I ld_A, char order_A, char op_A, char fill_A, T * X, I ld_X, char order_X, char op_X); \
    template void herk<T,I>(handle & h, I n, I k, T * A, I ld_A, char order_A, char op_A, T * C, I ld_C, char order_C, char fill_C); \
    template void hemv<T,I>(handle & h, I n, T * A, I ld_A, char order_A, char op_A, char fill_A, T * x, T * y); \
    template void gemv<T,I>(handle & h, I m, I n, T * A, I ld_A, char order_A, char op_A, T * x, T * y);

        #define INSTANTIATE_T(T) \
        INSTANTIATE_T_I(T, int32_t) \
        /* INSTANTIATE_T_I(T, int64_t) */

            // INSTANTIATE_T(float)
            INSTANTIATE_T(double)
            // INSTANTIATE_T(std::complex<float>)
            INSTANTIATE_T(std::complex<double>)

        #undef INSTANTIATE_T
    #undef INSTANTIATE_T_I

}
}
}
