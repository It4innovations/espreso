
#ifndef SRC_GPU_KERNELS_H_
#define SRC_GPU_KERNELS_H_

#include "gpu_management.h"
#include "math/primitives/vector_dense.h"
#include "math/primitives_new/vector_dense_view_new.h"

namespace espreso {
namespace gpu {
namespace kernels {

    template<typename T, typename I, typename A>
    void DCmap_scatter(mgm::queue & q, Vector_Dense<T*,I,A> & domain_vector_pointers, const Vector_Dense<I,I,A> & n_dofs_interfaces, const Vector_Dense<T,I,A> & cluster_vector, const Vector_Dense<I*,I,A> & D2Cs);

    template<typename T, typename I>
    void DCmap_scatter_new(mgm::queue & q, const VectorDenseView_new<T> & vec_cluster, MultiVectorDenseView_new<T,I> & vecs_subdomains, const MultiVectorDenseView_new<I,I> & D2C);

    template<typename T, typename I, typename A>
    void DCmap_gather(mgm::queue & q, const Vector_Dense<T*,I,A> & domain_vector_pointers, const Vector_Dense<I,I,A> & n_dofs_interfaces, Vector_Dense<T,I,A> & cluster_vector, const Vector_Dense<I*,I,A> & D2Cs);

    template<typename T, typename I>
    void DCmap_gather_new(mgm::queue & q, VectorDenseView_new<T> & vec_cluster, const MultiVectorDenseView_new<T,I> & vecs_subdomains, const MultiVectorDenseView_new<I,I> & D2C);

    template<typename T, typename I, typename Ao, typename Ai>
    void copy_matrix_triangle(mgm::queue & q, Matrix_Dense<T,I,Ao> & output, const Matrix_Dense<T,I,Ai> & input, char fill, char order);

}
}
}

#endif /* SRC_GPU_KERNELS_H_ */
