
namespace espreso {

    #define INSTANTIATE_T_I_M(T,I,M) \
    template struct SpBLAS<M,T,I>;

        #define INSTANTIATE_T_I(T,I) \
        INSTANTIATE_T_I_M(T,I,Matrix_CSR)

            #define INSTANTIATE_T(T) \
            INSTANTIATE_T_I(T,int32_t) \
            /* INSTANTIATE_T_I(T,int64_t) */

                INSTANTIATE_T(float)
                INSTANTIATE_T(double)
                // INSTANTIATE_T(std::complex<float>)
                INSTANTIATE_T(std::complex<double>)

            #undef INSTANTIATE_T
        #undef INSTANTIATE_T_I
    #undef INSTANTIATE_T_I_M

}
