
#if !defined(HAVE_CUDA) && !defined(HAVE_ROCM) && !defined(HAVE_ONEAPI)

#include "gpu_kernels.h"

namespace espreso {
namespace gpu {
namespace kernels {

    template<typename T, typename I, typename A>
    void DCmap_scatter(mgm::queue & q, Vector_Dense<T*,I,A> & domain_vector_pointers, const Vector_Dense<I,I,A> & n_dofs_interfaces, const Vector_Dense<T,I,A> & cluster_vector, const Vector_Dense<I*,I,A> & D2Cs) {}

    template<typename T, typename I, typename A>
    void DCmap_gather(mgm::queue & q, const Vector_Dense<T*,I,A> & domain_vector_pointers, const Vector_Dense<I,I,A> & n_dofs_interfaces, Vector_Dense<T,I,A> & cluster_vector, const Vector_Dense<I*,I,A> & D2Cs) {}

    template<typename T, typename I, typename Ao, typename Ai>
    void copy_matrix_triangle(mgm::queue & q, Matrix_Dense<T,I,Ao> & output, const Matrix_Dense<T,I,Ai> & input, char fill, char order) {}

}
}
}

#include "gpu/gpu_kernels.inst.hpp"

#endif
