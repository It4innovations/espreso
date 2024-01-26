
#ifdef HAVE_CUDA

#include "gpu/gpu_kernels.h"

#include "w.cuda.common.h"



namespace espreso {
namespace gpu {
namespace kernels {

    namespace
    {
        template<typename T, typename I>
        static __global__ void _do_DCmap_scatter(T ** domain_vectors, const I * n_dofs_interfaces, const T * cluster_vector, I const * const * D2Cs)
        {
            // one block per domain
        
            I d = blockIdx.x;
            I n_dofs_interface = n_dofs_interfaces[d];
            T * domain_vector = domain_vectors[d];
            const I * D2C = D2Cs[d];
        
            for(I dof = threadIdx.x; dof < n_dofs_interface; dof += blockDim.x)
            {
                domain_vector[dof] = cluster_vector[D2C[dof]];
            }
        }
        
        template<typename T, typename I>
        static __global__ void _do_DCmap_gather(T const * const * domain_vectors, const I * n_dofs_interfaces, T * cluster_vector, I const * const * D2Cs)
        {
            // one block per domain
        
            I d = blockIdx.x;
            I n_dofs_interface = n_dofs_interfaces[d];
            const T * domain_vector = domain_vectors[d];
            const I * D2C = D2Cs[d];
        
            for(I dof = threadIdx.x; dof < n_dofs_interface; dof += blockDim.x)
            {
                atomicAdd(&cluster_vector[D2C[dof]], domain_vector[dof]);
            }
        }
    }

    template<typename T, typename I>
    void DCmap_scatter(mgm::queue & q, Vector_Dense<T*,I,mgm::Ad> & domain_vector_pointers, const Vector_Dense<I,I,mgm::Ad> & n_dofs_interfaces, const Vector_Dense<T,I,mgm::Ad> & cluster_vector, const Vector_Dense<I*,I,mgm::Ad> & D2Cs)
    {
        I n_domains = domain_vector_pointers.size;
        _do_DCmap_scatter<T,I><<< n_domains, 256, 0, q.stream >>>(domain_vector_pointers.vals, n_dofs_interfaces.vals, cluster_vector.vals, D2Cs.vals);
        CHECK(cudaPeekAtLastError());
    }

    template<typename T, typename I>
    void DCmap_gather(mgm::queue & q, const Vector_Dense<T*,I,mgm::Ad> & domain_vector_pointers, const Vector_Dense<I,I,mgm::Ad> & n_dofs_interfaces, Vector_Dense<T,I,mgm::Ad> & cluster_vector, const Vector_Dense<I*,I,mgm::Ad> & D2Cs)
    {
        I n_domains = domain_vector_pointers.size;
        _do_DCmap_gather<T,I><<< n_domains, 256, 0, q.stream >>>(domain_vector_pointers.vals, n_dofs_interfaces.vals, cluster_vector.vals, D2Cs.vals);
        CHECK(cudaPeekAtLastError());
    }



    #define INSTANTIATE(T,I) \
    template void DCmap_scatter<T,I>(mgm::queue & q, Vector_Dense<T*,I,mgm::Ad> & domain_vector_pointers, const Vector_Dense<I,I,mgm::Ad> & n_dofs_interfaces, const Vector_Dense<T,I,mgm::Ad> & cluster_vector, const Vector_Dense<I*,I,mgm::Ad> & D2Cs); \
    template void DCmap_gather<T,I>(mgm::queue & q, const Vector_Dense<T*,I,mgm::Ad> & domain_vector_pointers, const Vector_Dense<I,I,mgm::Ad> & n_dofs_interfaces, Vector_Dense<T,I,mgm::Ad> & cluster_vector, const Vector_Dense<I*,I,mgm::Ad> & D2Cs);

    INSTANTIATE(float,  int32_t)
    INSTANTIATE(float,  int64_t)
    INSTANTIATE(double, int32_t)
    INSTANTIATE(double, int64_t)

    #undef INSTANTIATE
}
}
}




#endif
