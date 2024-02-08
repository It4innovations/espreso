
#ifdef HAVE_CUDA

#include "gpu/gpu_kernels.h"
#include "w.cuda.gpu_management.h"
#include "basis/utilities/utils.h"

#include <complex>



namespace espreso {
namespace gpu {
namespace kernels {

    namespace
    {
        template<typename T>
        static __device__ void complexAtomicAdd(std::complex<T> * dst, std::complex<T> val)
        {
            atomicAdd(&reinterpret_cast<T*>(dst)[0], reinterpret_cast<T*>(&val)[0]);
            atomicAdd(&reinterpret_cast<T*>(dst)[1], reinterpret_cast<T*>(&val)[1]);
        }

        template<typename T>
        static __device__ void myAtomicAdd(T * dst, T val) { atomicAdd(dst, val); }
        template<typename T>
        static __device__ void myAtomicAdd(std::complex<T> * dst, std::complex<T> val) { complexAtomicAdd(dst, val); }

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
            // launch with one block per domain
            #if defined(__CUDA_ARCH__) && __CUDA_ARCH__ >= 600
            #define MY_DOUBLE_ATOMIC_AVAILABLE true
            #else
            #define MY_DOUBLE_ATOMIC_AVAILABLE false
            #endif

            if constexpr(std::is_same_v<utils::remove_complex_t<T>,double> && !MY_DOUBLE_ATOMIC_AVAILABLE)
            {
                // must iterate sequentially through the domains
                I n_domains = blockDim.x;
                if(blockIdx.x > 0) return;

                for(I d = 0; d < n_domains; d++)
                {
                    I n_dofs_interface = n_dofs_interfaces[d];
                    const T * domain_vector = domain_vectors[d];
                    const I * D2C = D2Cs[d];
                
                    for(I dof = threadIdx.x; dof < n_dofs_interface; dof += blockDim.x)
                    {
                        cluster_vector[D2C[dof]] += domain_vector[dof];
                    }
                }
            }
            else
            {
                I d = blockIdx.x;
                I n_dofs_interface = n_dofs_interfaces[d];
                const T * domain_vector = domain_vectors[d];
                const I * D2C = D2Cs[d];
            
                for(I dof = threadIdx.x; dof < n_dofs_interface; dof += blockDim.x)
                {
                    myAtomicAdd(&cluster_vector[D2C[dof]], domain_vector[dof]);
                }
            }
        }
    }

    template<typename T, typename I>
    void DCmap_scatter(mgm::queue & q, Vector_Dense<T*,I,mgm::Ad> & domain_vector_pointers, const Vector_Dense<I,I,mgm::Ad> & n_dofs_interfaces, const Vector_Dense<T,I,mgm::Ad> & cluster_vector, const Vector_Dense<I*,I,mgm::Ad> & D2Cs)
    {
        I n_domains = domain_vector_pointers.size;
        _do_DCmap_scatter<T,I><<< n_domains, 256, 0, q->stream >>>(domain_vector_pointers.vals, n_dofs_interfaces.vals, cluster_vector.vals, D2Cs.vals);
        CHECK(cudaPeekAtLastError());
    }

    template<typename T, typename I>
    void DCmap_gather(mgm::queue & q, const Vector_Dense<T*,I,mgm::Ad> & domain_vector_pointers, const Vector_Dense<I,I,mgm::Ad> & n_dofs_interfaces, Vector_Dense<T,I,mgm::Ad> & cluster_vector, const Vector_Dense<I*,I,mgm::Ad> & D2Cs)
    {
        I n_domains = domain_vector_pointers.size;
        _do_DCmap_gather<T,I><<< n_domains, 256, 0, q->stream >>>(domain_vector_pointers.vals, n_dofs_interfaces.vals, cluster_vector.vals, D2Cs.vals);
        CHECK(cudaPeekAtLastError());
    }



    #define INSTANTIATE_T_I(T,I) \
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

}
}
}

#endif
