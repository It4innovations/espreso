
#ifdef HAVE_CUDA
#ifdef ESPRESO_USE_WRAPPER_GPU_CUDA

#include "gpu/gpu_kernels.h"
#include "common_cuda_mgm.h"
#include "common_internal.cuh"
#include "basis/utilities/utils.h"

#include <complex>



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
        
            for(I dof = threadIdx.x; dof < n_dofs_interface; dof += blockDim.x) {
                domain_vector[dof] = cluster_vector[D2C[dof]];
            }
        }

        template<typename T, typename I>
        static __global__ void _do_DCmap_scatter_new(const T * vec_cluster, const I * vec_subdomains_offsets, T * vec_subdomains_data, const I * D2C_offsets, const I * D2C_data)
        {
            // one block per domain
        
            I di = blockIdx.x;
            I start = vec_subdomains_offsets[di];
            I end = vec_subdomains_offsets[di+1];
            I size = end - start;
            T * vec = vec_subdomains_data + start;
            const I * D2C = D2C_data + D2C_offsets[di];
        
            for(I i = threadIdx.x; i < size; i += blockDim.x) {
                vec[i] = vec_cluster[D2C[i]];
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

            if constexpr(std::is_same_v<utils::remove_complex_t<T>,double> && !MY_DOUBLE_ATOMIC_AVAILABLE) {
                // must iterate sequentially through the domains
                I n_domains = blockDim.x;
                if(blockIdx.x > 0) return;

                for(I d = 0; d < n_domains; d++) {
                    I n_dofs_interface = n_dofs_interfaces[d];
                    const T * domain_vector = domain_vectors[d];
                    const I * D2C = D2Cs[d];
                
                    for(I dof = threadIdx.x; dof < n_dofs_interface; dof += blockDim.x) {
                        myGenericAdd(cluster_vector + D2C[dof], domain_vector[dof]);
                    }
                }
            }
            else {
                I d = blockIdx.x;
                I n_dofs_interface = n_dofs_interfaces[d];
                const T * domain_vector = domain_vectors[d];
                const I * D2C = D2Cs[d];
            
                for(I dof = threadIdx.x; dof < n_dofs_interface; dof += blockDim.x) {
                    myAtomicAdd(&cluster_vector[D2C[dof]], domain_vector[dof]);
                }
            }
        }

        template<typename T, typename I>
        static __global__ void _do_DCmap_gather_new(T * vec_cluster, const I * vec_subdomains_offsets, const T * vec_subdomains_data, const I * D2C_offsets, const I * D2C_data)
        {
            // launch with one block per domain

            #if defined(__CUDA_ARCH__) && __CUDA_ARCH__ >= 600
            #define MY_DOUBLE_ATOMIC_AVAILABLE true
            #else
            #define MY_DOUBLE_ATOMIC_AVAILABLE false
            #endif

            if constexpr(std::is_same_v<utils::remove_complex_t<T>,double> && !MY_DOUBLE_ATOMIC_AVAILABLE) {
                // iterate sequentially through the domains
                I n_domains = blockDim.x;
                if(blockIdx.x > 0) return;

                for(I di = 0; di < n_domains; di++) {
                    I start = vec_subdomains_offsets[di];
                    I end = vec_subdomains_offsets[di+1];
                    I size = end - start;
                    const T * vec = vec_subdomains_data + start;
                    const I * D2C = D2C_data + D2C_offsets[di];
                
                    for(I i = threadIdx.x; i < size; i += blockDim.x) {
                        myGenericAdd(vec_cluster + D2C[i], vec[i]);
                    }
                }
            }
            else {
                int di = blockIdx.x;
                I start = vec_subdomains_offsets[di];
                I end = vec_subdomains_offsets[di+1];
                I size = end - start;
                const T * vec = vec_subdomains_data + start;
                const I * D2C = D2C_data + D2C_offsets[di];
            
                for(I i = threadIdx.x; i < size; i += blockDim.x) {
                    myAtomicAdd(&vec_cluster[D2C[i]], vec[i]);
                }
            }
        }

        template<typename T, typename I>
        static __global__ void _copy_matrix_triangle(T * output, I output_ld, T const * input, I input_ld, I n, char fill)
        {
            // assuming row-major matrix
            // launch with (n-1)/2+1 blocks (one block per two rows)
            // each block handles rows blockIdx.x (for stage=0) and n-blockIdx.x-1 (for stage=1).
            // in stage=0, warps are normally sequential next to each other
            // in stage=1, the order of warps is reversed for better load balance between warps

            I warpidx = threadIdx.x / warpSize;
            I warpcount = blockDim.x / warpSize;
            I idxinwarp = threadIdx.x % warpSize;
            #pragma unroll
            for(int stage = 0; stage <= 1; stage++) {
                if(stage == 1 && blockIdx.x == n/2) return;
                I r = (stage == 0) ? (blockIdx.x) : (n - blockIdx.x - 1);
                I chunk_idx = (stage == 0) ? (warpidx) : (warpcount - warpidx - 1);
                I myidx = chunk_idx * warpSize + idxinwarp;
                const T * row_in = input + r * input_ld;
                T * row_out = output + r * output_ld;
                if(fill == 'L') {
                    for(I c = myidx; c <= r; c += blockDim.x) row_out[c] = row_in[c];
                }
                if(fill == 'U') {
                    I skip_leaps = r / blockDim.x;
                    for(I c = skip_leaps * blockDim.x + myidx; c < n; c += blockDim.x) if(c >= r) row_out[c] = row_in[c];
                }
            }
        }
    }

    template<typename T, typename I, typename A>
    void DCmap_scatter(mgm::queue & q, Vector_Dense<T*,I,A> & domain_vector_pointers, const Vector_Dense<I,I,A> & n_dofs_interfaces, const Vector_Dense<T,I,A> & cluster_vector, const Vector_Dense<I*,I,A> & D2Cs)
    {
        static_assert(A::is_data_device_accessible, "data has to be device accessible");
        I n_domains = domain_vector_pointers.size;
        _do_DCmap_scatter<T,I><<< n_domains, 256, 0, q->stream >>>(domain_vector_pointers.vals, n_dofs_interfaces.vals, cluster_vector.vals, D2Cs.vals);
        CHECK(cudaPeekAtLastError());
    }

    template<typename T, typename I>
    void DCmap_scatter_new(mgm::queue & q, const VectorDenseView_new<T> & vec_cluster, MultiVectorDenseView_new<T,I> & vecs_subdomains, const MultiVectorDenseView_new<I,I> & D2C)
    {
        if(!vec_cluster.ator->is_data_accessible_gpu()) eslog::error("wrong allocator\n");
        if(!vecs_subdomains.ator->is_data_accessible_gpu()) eslog::error("wrong allocator\n");
        if(!D2C.ator->is_data_accessible_gpu()) eslog::error("wrong allocator\n");
        I n_domains = vecs_subdomains.num_vectors;
        _do_DCmap_scatter_new<T,I><<< n_domains, 256, 0, q->stream >>>(vec_cluster.vals, vecs_subdomains.offsets, vecs_subdomains.vals, D2C.offsets, D2C.vals);
        CHECK(cudaPeekAtLastError());
    }

    template<typename T, typename I, typename A>
    void DCmap_gather(mgm::queue & q, const Vector_Dense<T*,I,A> & domain_vector_pointers, const Vector_Dense<I,I,A> & n_dofs_interfaces, Vector_Dense<T,I,A> & cluster_vector, const Vector_Dense<I*,I,A> & D2Cs)
    {
        static_assert(A::is_data_device_accessible, "data has to be device accessible");
        I n_domains = domain_vector_pointers.size;
        _do_DCmap_gather<T,I><<< n_domains, 256, 0, q->stream >>>(domain_vector_pointers.vals, n_dofs_interfaces.vals, cluster_vector.vals, D2Cs.vals);
        CHECK(cudaPeekAtLastError());
    }

    template<typename T, typename I>
    void DCmap_gather_new(mgm::queue & q, VectorDenseView_new<T> & vec_cluster, const MultiVectorDenseView_new<T,I> & vecs_subdomains, const MultiVectorDenseView_new<I,I> & D2C)
    {
        if(!vec_cluster.ator->is_data_accessible_gpu()) eslog::error("wrong allocator\n");
        if(!vecs_subdomains.ator->is_data_accessible_gpu()) eslog::error("wrong allocator\n");
        if(!D2C.ator->is_data_accessible_gpu()) eslog::error("wrong allocator\n");
        I n_domains = vecs_subdomains.num_vectors;
        _do_DCmap_gather_new<T,I><<< n_domains, 256, 0, q->stream >>>(vec_cluster.vals, vecs_subdomains.offsets, vecs_subdomains.vals, D2C.offsets, D2C.vals);
        CHECK(cudaPeekAtLastError());
    }

    template<typename T, typename I, typename Ao, typename Ai>
    void copy_matrix_triangle(mgm::queue & q, Matrix_Dense<T,I,Ao> & output, const Matrix_Dense<T,I,Ai> & input, char fill, char order)
    {
        static_assert(Ao::is_data_device_accessible, "matrix data has to be device accessible");
        static_assert(Ai::is_data_device_accessible, "matrix data has to be device accessible");
        if(output.nrows != input.nrows || output.ncols != input.ncols || input.nrows != input.ncols) eslog::error("matrix dimensions do not match\n");

        if(order == 'R') {
            int bpg = (input.nrows - 1) / 2 + 1;
            _copy_matrix_triangle<<< bpg, 256, 0, q->stream >>>(output.vals, output.get_ld(), input.vals, input.get_ld(), input.nrows, fill);
            CHECK(cudaPeekAtLastError());
        }
        else if(order == 'C') {
            char fill_compl = mgm::fill_change(fill);
            copy_matrix_triangle(q, output, input, fill_compl, 'R');
        }
        else {
            eslog::error("invalid order %c\n", order);
        }
    }

}
}
}

#include "gpu/gpu_kernels.inst.hpp"

#endif
#endif
