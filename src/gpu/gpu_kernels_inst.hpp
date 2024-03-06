
#include "basis/utilities/cbmb_allocator.h"

namespace espreso {
namespace gpu {
namespace kernels {

    #define INSTANTIATE_T_I_AO_AI(T,I,AO,AI) \
    template void copy_matrix_triangle<T,I,AO,AI>(mgm::queue & q, Matrix_Dense<T,I,AO> & output, const Matrix_Dense<T,I,AI> & input, char fill, char order);

        #define INSTANTIATE_T_I_AO(T,I,AO) \
        INSTANTIATE_T_I_AO_AI(T,I,AO,mgm::Ad) \
        INSTANTIATE_T_I_AO_AI(T,I,AO,cbmba_d)

            #define INSTANTIATE_T_I(T,I) \
            INSTANTIATE_T_I_AO(T,I,mgm::Ad) \
            INSTANTIATE_T_I_AO(T,I,cbmba_d) \
            template void DCmap_scatter<T,I>(mgm::queue & q, Vector_Dense<T*,I,mgm::Ad> & domain_vector_pointers, const Vector_Dense<I,I,mgm::Ad> & n_dofs_interfaces, const Vector_Dense<T,I,mgm::Ad> & cluster_vector, const Vector_Dense<I*,I,mgm::Ad> & D2Cs); \
            template void DCmap_gather<T,I>(mgm::queue & q, const Vector_Dense<T*,I,mgm::Ad> & domain_vector_pointers, const Vector_Dense<I,I,mgm::Ad> & n_dofs_interfaces, Vector_Dense<T,I,mgm::Ad> & cluster_vector, const Vector_Dense<I*,I,mgm::Ad> & D2Cs);

                #define INSTANTIATE_T(T) \
                INSTANTIATE_T_I(T, int32_t) \
                /* INSTANTIATE_T_I(T, int64_t) */

                    // INSTANTIATE_T(float)
                    INSTANTIATE_T(double)
                    // INSTANTIATE_T(std::complex<float>)
                    // INSTANTIATE_T(std::complex<double>)

                #undef INSTANTIATE_T
            #undef INSTANTIATE_T_I
        #undef INSTANTIATE_T_I_AO
    #undef INSTANTIATE_T_I_AO_AI

}
}
}
