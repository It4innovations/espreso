
namespace espreso {

    #define INSTANTIATE_T_I(T,I) \
    template struct DenseSolver<T,I>;

        #define INSTANTIATE_T(T) \
        INSTANTIATE_T_I(T,int32_t); \
        /* INSTANTIATE_T_I(T,int64_t); */

            // INSTANTIATE_T(float)
            INSTANTIATE_T(double)
            // INSTANTIATE_T(std::complex<float>)
            INSTANTIATE_T(std::complex<double>)

        #undef INSTANTIATE_T
    #undef INSTANTIATE_T_I

}
