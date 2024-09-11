
#include "math.spblas.h"

namespace espreso {

    #define INSTANTIATE_T_I_M(T,I,M) \
    template struct SpBLAS<M,T,I>;

        #define INSTANTIATE_T_I(T,I) \
        INSTANTIATE_T_I_M(T,I,Matrix_CSR)

            #define INSTANTIATE_T(T) \
            INSTANTIATE_T_I(T,int32_t) \
            INSTANTIATE_T_I(T,int64_t)

                INSTANTIATE_T(float)
                INSTANTIATE_T(double)
                // INSTANTIATE_T(std::complex<float>)
                INSTANTIATE_T(std::complex<double>)

            #undef INSTANTIATE_T
        #undef INSTANTIATE_T_I
    #undef INSTANTIATE_T_I_M


namespace math {
namespace spblas {

#define INSTANTIATE_SUBMATRIX(T, I) \
template void submatrix<T, I>(const Matrix_CSR<T, I>&, Matrix_Dense<T, I>&, I, I, I, I, bool, bool, bool); \
template void submatrix<T, I>(const Matrix_CSR<T, I>&, Matrix_CSR<T, I>&  , I, I, I, I, bool, bool, bool);

INSTANTIATE_SUBMATRIX(double, int)
INSTANTIATE_SUBMATRIX(std::complex<double>, int)

#undef INSTANTIATE_T_I

}
}
}
