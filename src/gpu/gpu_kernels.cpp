
#if !defined(HAVE_CUDA) && !defined(HAVE_ROCM)

#include "gpu_kernels.h"

namespace espreso {
namespace gpu {
namespace kernels {

    template<typename T, typename I>
    void DCmap_scatter(mgm::queue & q, Vector_Dense<T*,I,mgm::Ad> & domain_vector_pointers, const Vector_Dense<I,I,mgm::Ad> & n_dofs_interfaces, const Vector_Dense<T,I,mgm::Ad> & cluster_vector, const Vector_Dense<I*,I,mgm::Ad> & D2Cs) {}

    template<typename T, typename I>
    void DCmap_gather(mgm::queue & q, const Vector_Dense<T*,I,mgm::Ad> & domain_vector_pointers, const Vector_Dense<I,I,mgm::Ad> & n_dofs_interfaces, Vector_Dense<T,I,mgm::Ad> & cluster_vector, const Vector_Dense<I*,I,mgm::Ad> & D2Cs) {}



    #define INSTANTIATE(T,I) \
    template void DCmap_scatter<T,I>(mgm::queue & q, Vector_Dense<T*,I,mgm::Ad> & domain_vector_pointers, const Vector_Dense<I,I,mgm::Ad> & n_dofs_interfaces, const Vector_Dense<T,I,mgm::Ad> & cluster_vector, const Vector_Dense<I*,I,mgm::Ad> & D2Cs); \
    template void DCmap_gather<T,I>(mgm::queue & q, const Vector_Dense<T*,I,mgm::Ad> & domain_vector_pointers, const Vector_Dense<I,I,mgm::Ad> & n_dofs_interfaces, Vector_Dense<T,I,mgm::Ad> & cluster_vector, const Vector_Dense<I*,I,mgm::Ad> & D2Cs);
        // INSTANTIATE(float,                int32_t)
        INSTANTIATE(double,               int32_t)
        // INSTANTIATE(std::complex<float >, int32_t)
        // INSTANTIATE(std::complex<double>, int32_t)
        // INSTANTIATE(float,                int64_t)
        // INSTANTIATE(double,               int64_t)
        // INSTANTIATE(std::complex<float >, int64_t)
        // INSTANTIATE(std::complex<double>, int64_t)
    #undef INSTANTIATE
}
}
}

#endif
