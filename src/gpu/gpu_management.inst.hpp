
#include "basis/utilities/cbmb_allocator.h"
#include "basis/utilities/arena_allocator.h"

namespace espreso {
namespace gpu {
namespace mgm {

    #define INSTANTIATE_T_I_Ao_Ai(T,I,Ao,Ai) \
    template void copy_submit<T,I,Ao,Ai>(queue & q, Vector_Dense<T,I,Ao> & output, const Vector_Dense<T,I,Ai> & input); \
    template void copy_submit<T,I,Ao,Ai>(queue & q, Matrix_Dense<T,I,Ao> & output, const Matrix_Dense<T,I,Ai> & input); \
    template void copy_submit<T,I,Ao,Ai>(queue & q, Matrix_CSR<T,I,Ao>   & output, const Matrix_CSR<T,I,Ai>   & input, bool copy_pattern, bool copy_vals); \

        #define INSTANTIATE_T_I_A1_A2(T,I,A1,A2) \
        INSTANTIATE_T_I_Ao_Ai(T, I, A1, A2) \
        INSTANTIATE_T_I_Ao_Ai(T, I, A2, A1) \
        /* INSTANTIATE_T_I_Ao_Ai(T, I, A1, A1) */ \
        /* INSTANTIATE_T_I_Ao_Ai(T, I, A2, A2) */

            #define INSTANTIATE_T_I(T,I) \
            INSTANTIATE_T_I_A1_A2(T, I, mgm::Ah,       mgm::Ad) \
            INSTANTIATE_T_I_A1_A2(T, I, mgm::Ah,       arena_d) \
            INSTANTIATE_T_I_A1_A2(T, I, mgm::Ah,       cbmba_d) \
            INSTANTIATE_T_I_A1_A2(T, I, cpu_allocator, arena_d) \
            INSTANTIATE_T_I_A1_A2(T, I, cpu_allocator, cbmba_d)

                #define INSTANTIATE_T(T) \
                template void copy_submit<T>(queue & q, T * dst, T const * src, size_t num_elements); \
                INSTANTIATE_T_I(T, int32_t) \
                /* INSTANTIATE_T_I(T, int64_t) */

                    // INSTANTIATE_T(float)
                    INSTANTIATE_T(double)
                    // INSTANTIATE_T(std::complex<float>)
                    // INSTANTIATE_T(std::complex<double>)
                    // INSTANTIATE_T(float*)
                    INSTANTIATE_T(double*)
                    // INSTANTIATE_T(std::complex<float>*)
                    // INSTANTIATE_T(std::complex<double>*)
                    INSTANTIATE_T(int32_t)
                    INSTANTIATE_T(int32_t*)
                    // INSTANTIATE_T(int64_t)
                    // INSTANTIATE_T(int64_t*)

                #undef INSTANTIATE_T
            #undef INSTANTIATE_T_I
        #undef INSTANTIATE_T_I_A1_A2
    #undef INSTANTIATE_T_I_Ao_Ai

}
}
}
