
#include "gpu_kernels.h"

namespace espreso {
namespace gpu {
namespace kernels {

#ifndef HAVE_CUDA

template<>
void DCmap_scatter(mgm::queue & q, Vector_Dense<double*,int,mgm::Ad> & domain_vector_pointers, const Vector_Dense<int,int,mgm::Ad> & n_dofs_interfaces, const Vector_Dense<double,int,mgm::Ad> & cluster_vector, const Vector_Dense<int*,int,mgm::Ad> & D2Cs) {}

template<>
void DCmap_gather(mgm::queue & q, const Vector_Dense<double*,int,mgm::Ad> & domain_vector_pointers, const Vector_Dense<int,int,mgm::Ad> & n_dofs_interfaces, Vector_Dense<double,int,mgm::Ad> & cluster_vector, const Vector_Dense<int*,int,mgm::Ad> & D2Cs) {}

#endif

}
}
}
