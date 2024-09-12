
#include "gpu/gpu_management.h"

namespace espreso {

    #define INSTANTIATE_T_I_A(T,I,A) \
    template void SchurComplementSolver<T, I>::factorizeNumericAndGetSc<A>(Matrix_Dense<T,I,A> & sc, char uplo, T alpha);

        #define INSTANTIATE_T_I(T,I) \
        INSTANTIATE_T_I_A(T,I,cpu_allocator) \
        INSTANTIATE_T_I_A(T,I,gpu::mgm::Ah) \
        template class SchurComplementSolver<T,I>;

            #define INSTANTIATE_T(T) \
            INSTANTIATE_T_I(T,int32_t); \
            /* INSTANTIATE_T_I(T,int64_t); */

                // INSTANTIATE_T(float)
                INSTANTIATE_T(double)
                // INSTANTIATE_T(std::complex<float>)
                // INSTANTIATE_T(std::complex<double>)

            #undef INSTANTIATE_T
        #undef INSTANTIATE_T_I
    #undef INSTANTIATE_T_I_A

}
