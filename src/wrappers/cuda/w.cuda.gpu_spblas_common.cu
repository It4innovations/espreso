
#ifdef HAVE_CUDA
#ifdef ESPRESO_USE_WRAPPER_GPU_CUDA

#include "w.cuda.gpu_spblas_common.h"
#include "w.cuda.gpu_management.h"

#include <type_traits>
#include <algorithm>
#include <complex>
#include <cub/cub.cuh>



namespace espreso {
namespace gpu {
namespace spblas {

    namespace {
        template<typename I>
        I get_most_significant_bit(I val)
        {
            static_assert(std::is_integral_v<I> && std::is_unsigned_v<I>, "wrong type");
            I msb = 0;
            while(val != 0) {
                val >>= 1;
                msb++;
            }
            return msb;
        }

        template<typename T> __device__ constexpr bool is_complex();
        template<> __device__ constexpr bool is_complex<int32_t>() { return false; }
        // template<> __device__ constexpr bool is_complex<int64_t>() { return false; }
        // template<> __device__ constexpr bool is_complex<float>() { return false; }
        template<> __device__ constexpr bool is_complex<double>() { return false; }
        // template<> __device__ constexpr bool is_complex<std::complex<float>>() { return true; }
        // template<> __device__ constexpr bool is_complex<std::complex<double>>() { return true; }

        template<typename T>
        static __global__ void _init_linear(T * output, size_t count)
        {
            size_t idx = blockIdx.x * blockDim.x + threadIdx.x;
            size_t stride = blockDim.x * gridDim.x;
            for(size_t i = idx; i < count; i += stride) output[i] = (T)i;
        }

        template<typename T, typename I, bool conj = false>
        static __global__ void _permute_array(T * output, T const * input, I const * perm, size_t count)
        {
            size_t idx = blockIdx.x * blockDim.x + threadIdx.x;
            size_t stride = blockDim.x * gridDim.x;
            for(size_t i = idx; i < count; i += stride) {
                T & out = output[i];
                const T & val = input[perm[i]];
                if constexpr(conj && is_complex<T>()) {
                    reinterpret_cast<T*>(&out)[0] =  reinterpret_cast<T*>(&val)[0];
                    reinterpret_cast<T*>(&out)[1] = -reinterpret_cast<T*>(&val)[1];
                }
                else {
                    out = val;
                }
            }
        }

        template<typename I>
        static __global__ void _csr_to_ijv_rowidxs(I * out_rowidxs, I const * in_rowptrs)
        {
            I r = blockIdx.x;
            I start = in_rowptrs[r];
            I end = in_rowptrs[r+1];
            for(I i = start + threadIdx.x; i < end; i += blockDim.x) out_rowidxs[i] = r;
        }

        template<typename I>
        static __global__ void _ijv_rowidxs_to_csr_rowptrs(I * csr_rowptrs, I * ijv_rowidxs_sorted, I nrows, I nnz)
        {
            I idx = blockIdx.x * blockDim.x + threadIdx.x;
            I stride = blockDim.x * gridDim.x;
            for(I i = idx; i < nnz-1; i += stride) {
                I curr_row = ijv_rowidxs_sorted[i];
                I next_row = ijv_rowidxs_sorted[i+1];
                for(I r = curr_row; r < next_row; r++) csr_rowptrs[r+1] = i+1;
            }
            if(idx == stride-1) {
                I lastrow = ijv_rowidxs_sorted[nnz-1];
                for(I r = lastrow; r < nrows; r++) csr_rowptrs[r+1] = nnz;

                I firstrow = ijv_rowidxs_sorted[0];
                for(I r = 0; r <= firstrow; r++) csr_rowptrs[r] = 0;
            }
        }
    }

    template<typename I>
    void my_csr_transpose_buffersize(cudaStream_t & stream, I input_nrows, I input_ncols, I nnz, size_t & buffersize)
    {
        I end_bit = get_most_significant_bit((uint64_t)input_ncols);
        size_t bfs_map_and_linear = nnz * sizeof(I);
        size_t bfs_sorted_colidxs = nnz * sizeof(I);
        size_t bfs_sort;
        CHECK(cub::DeviceRadixSort::SortPairs(nullptr, bfs_sort, (I*)nullptr, (I*)nullptr, (I*)nullptr, (I*)nullptr, nnz, 0, end_bit, stream));
        CHECK(cudaStreamSynchronize(stream));
        buffersize = bfs_map_and_linear + bfs_sorted_colidxs + bfs_sort;
    }
    
    template<typename I>
    void my_csr_transpose_preprocess(cudaStream_t & stream, I input_nrows, I input_ncols, I nnz, I * input_rowptrs, I * input_colidxs, I * output_rowptrs, I * output_colidxs, size_t buffersize, void * buffer)
    {
        I output_nrows = input_ncols;
        I end_bit = get_most_significant_bit((uint64_t)input_ncols);
        I * map = (I*)buffer;
        buffer = (char*)buffer + nnz * sizeof(I);   buffersize -= nnz * sizeof(I);
        I * colidxs_sorted = (I*)buffer;
        buffer = (char*)buffer + nnz * sizeof(I);   buffersize -= nnz * sizeof(I);
        I * ijv_rowidxs = colidxs_sorted; // just two unrelated temporary buffers sharing the same memory
        I * linear = map; // also sharing memory
        _init_linear<<< 16, 256, 0, stream >>>(linear, nnz);
        CHECK(cudaPeekAtLastError());
        CHECK(cub::DeviceRadixSort::SortPairs(buffer, buffersize, input_colidxs, colidxs_sorted, linear, map, nnz, 0, end_bit, stream));
        _ijv_rowidxs_to_csr_rowptrs<<< 16, 256, 0, stream >>>(output_rowptrs, colidxs_sorted, output_nrows, nnz);
        CHECK(cudaPeekAtLastError());
        _csr_to_ijv_rowidxs<<< input_nrows, 64, 0, stream >>>(ijv_rowidxs, input_rowptrs);
        CHECK(cudaPeekAtLastError());
        _permute_array<<< 16, 256, 0, stream >>>(output_colidxs, ijv_rowidxs, map, nnz);
        CHECK(cudaPeekAtLastError());
    }
    
    template<typename T, typename I>
    void my_csr_transpose_compute(cudaStream_t & stream, I nnz, T * input_vals, T * output_vals, bool conjugate, void * buffer)
    {
        I * map = (I*)buffer;
        if(conjugate) _permute_array<T,I,true> <<< 16, 256, 0, stream >>>(output_vals, input_vals, map, nnz);
        else          _permute_array<T,I,false><<< 16, 256, 0, stream >>>(output_vals, input_vals, map, nnz);
        CHECK(cudaPeekAtLastError());
    }



    #define INSTANTIATE_T_I(T,I) \
    template void my_csr_transpose_compute<T,I>(cudaStream_t & stream, I nnz, T * input_vals, T * output_vals, bool conjugate, void * buffer);

        #define INSTANTIATE_I(I) \
        /* INSTANTIATE_T_I(float,I) */ \
        INSTANTIATE_T_I(double,I) \
        /* INSTANTIATE_T_I(std::complex<float>,I) */ \
        /* INSTANTIATE_T_I(std::complex<double>,I) */ \
        template void my_csr_transpose_buffersize<I>(cudaStream_t & stream, I input_nrows, I input_ncols, I nnz, size_t & buffersize); \
        template void my_csr_transpose_preprocess<I>(cudaStream_t & stream, I input_nrows, I input_ncols, I nnz, I * input_rowptrs, I * input_colidxs, I * output_rowptrs, I * output_colidxs, size_t buffersize, void * buffer);

            INSTANTIATE_I(int32_t)
            // INSTANTIATE_I(int64_t)

        #undef INSTANTIATE_I
    #undef INSTANTIATE_T_I

}
}
}

#endif
#endif
