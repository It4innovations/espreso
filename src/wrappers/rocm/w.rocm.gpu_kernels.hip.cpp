
#ifdef HAVE_ROCM

#include "gpu/gpu_kernels.h"
#include "w.rocm.gpu_management.h"

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
        __global__ void _do_DCmap_scatter(T ** domain_vectors, const I * n_dofs_interfaces, const T * cluster_vector, I const * const * D2Cs)
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
        __global__ void _do_DCmap_gather(T const * const * domain_vectors, const I * n_dofs_interfaces, T * cluster_vector, I const * const * D2Cs)
        {
            // one block per domain

            I d = blockIdx.x;
            I n_dofs_interface = n_dofs_interfaces[d];
            const T * domain_vector = domain_vectors[d];
            const I * D2C = D2Cs[d];

            for(I dof = threadIdx.x; dof < n_dofs_interface; dof += blockDim.x)
            {
                myAtomicAdd(&cluster_vector[D2C[dof]], domain_vector[dof]);
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
            for(int stage = 0; stage <= 1; stage++)
            {
                if(stage == 1 && blockIdx.x == n/2) return;
                I r = (stage == 0) ? (blockIdx.x) : (n - blockIdx.x - 1);
                I chunk_idx = (stage == 0) ? (warpidx) : (warpcount - warpidx - 1);
                I myidx = chunk_idx * warpSize + idxinwarp;
                const T * row_in = input + r * input_ld;
                T * row_out = output + r * output_ld;
                if(fill == 'L')
                {
                    for(I c = myidx; c <= r; c += blockDim.x) row_out[c] = row_in[c];
                }
                if(fill == 'U')
                {
                    I skip_leaps = r / blockDim.x;
                    for(I c = skip_leaps * blockDim.x + myidx; c < n; c += blockDim.x) if(c >= r) row_out[c] = row_in[c];
                }
            }
        }
    }

    template<typename T, typename I>
    void DCmap_scatter(mgm::queue & q, Vector_Dense<T*,I,mgm::Ad> & domain_vector_pointers, const Vector_Dense<I,I,mgm::Ad> & n_dofs_interfaces, const Vector_Dense<T,I,mgm::Ad> & cluster_vector, const Vector_Dense<I*,I,mgm::Ad> & D2Cs)
    {
        I n_domains = domain_vector_pointers.size;
        _do_DCmap_scatter<T,I><<< n_domains, 256, 0, q->stream >>>(domain_vector_pointers.vals, n_dofs_interfaces.vals, cluster_vector.vals, D2Cs.vals);
        CHECK(hipPeekAtLastError());
    }

    template<typename T, typename I>
    void DCmap_gather(mgm::queue & q, const Vector_Dense<T*,I,mgm::Ad> & domain_vector_pointers, const Vector_Dense<I,I,mgm::Ad> & n_dofs_interfaces, Vector_Dense<T,I,mgm::Ad> & cluster_vector, const Vector_Dense<I*,I,mgm::Ad> & D2Cs)
    {
        I n_domains = domain_vector_pointers.size;
        _do_DCmap_gather<T,I><<< n_domains, 256, 0, q->stream >>>(domain_vector_pointers.vals, n_dofs_interfaces.vals, cluster_vector.vals, D2Cs.vals);
        CHECK(hipPeekAtLastError());
    }

    template<typename T, typename I, typename Ao, typename Ai>
    void copy_matrix_triangle(mgm::queue & q, Matrix_Dense<T,I,Ao> & output, const Matrix_Dense<T,I,Ai> & input, char fill, char order)
    {
        static_assert(Ao::is_data_device_accessible, "matrix data has to be device accessible");
        static_assert(Ai::is_data_device_accessible, "matrix data has to be device accessible");
        if(output.nrows != input.nrows || output.ncols != input.ncols || input.nrows != input.ncols) eslog::error("matrix dimensions do not match\n");

        if(order == 'R')
        {
            int bpg = (input.nrows - 1) / 2 + 1;
            _copy_matrix_triangle<<< bpg, 256, 0, q->stream >>>(output.vals, output.get_ld(), input.vals, input.get_ld(), input.nrows, fill);
            CHECK(hipPeekAtLastError());
        }
        else if(order == 'C')
        {
            char fill_compl = mgm::fill_change(fill);
            copy_matrix_triangle(q, output, input, fill_compl, 'R');
        }
        else eslog::error("invalid order %c\n", order);
    }

}
}
}

#include "gpu/gpu_kernels.inst.hpp"

#endif
