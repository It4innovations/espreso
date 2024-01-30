
#ifndef SRC_GPU_KERNELS_H_
#define SRC_GPU_KERNELS_H_

#include "gpu_management.h"
#include "math/primitives/vector_dense.h"

namespace espreso {
namespace gpu {
namespace kernels {

    template<typename T, typename I>
    void DCmap_scatter(mgm::queue & q, Vector_Dense<T*,I,mgm::Ad> & domain_vector_pointers, const Vector_Dense<I,I,mgm::Ad> & n_dofs_interfaces, const Vector_Dense<T,I,mgm::Ad> & cluster_vector, const Vector_Dense<I*,I,mgm::Ad> & D2Cs);

    template<typename T, typename I>
    void DCmap_gather(mgm::queue & q, const Vector_Dense<T*,I,mgm::Ad> & domain_vector_pointers, const Vector_Dense<I,I,mgm::Ad> & n_dofs_interfaces, Vector_Dense<T,I,mgm::Ad> & cluster_vector, const Vector_Dense<I*,I,mgm::Ad> & D2Cs);

}
}
}

#endif /* SRC_GPU_KERNELS_H_ */
