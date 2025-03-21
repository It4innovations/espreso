
#ifdef HAVE_ROCM
#ifdef ESPRESO_USE_WRAPPER_GPU_ROCM

#include "gpu/gpu_management.h"
#include "w.rocm.gpu_management.h"
#include "basis/utilities/cbmb_allocator.h"

#include <omp.h>
#include <complex>



namespace espreso {
namespace gpu {
namespace mgm {

    gpu_wrapper_impl get_implementation()
    {
        return gpu_wrapper_impl::ROCM;
    }

    device get_device_by_mpi(int mpi_rank, int mpi_size)
    {
#ifndef ESPRESO_RANK_TO_GPU_MAP
#error "Undefined macro ESPRESO_RANK_TO_GPU_MAP. It should be defined in some wscript"
#endif
        static constexpr int rank_gpu_map[] = {ESPRESO_RANK_TO_GPU_MAP};
        static constexpr int n_gpus = sizeof(rank_gpu_map) / sizeof(*rank_gpu_map);

        device d = std::make_shared<_device>();
        if(mpi_size == 1) {
            d->gpu_idx = 0;
        }
        else if(mpi_size % n_gpus == 0) {
            // assuming rank0=gpu0, rank1=gpu1 etc..., and there are as many ranks on a node as gpus
            int local_node_rank = mpi_rank % n_gpus;
            d->gpu_idx = rank_gpu_map[local_node_rank];
            printf("rank %d localrank %d assigned gpu %d. btw ngpus=%d\n", mpi_rank, local_node_rank, d->gpu_idx, n_gpus); fflush(stdout);
        }
        else {
            eslog::error("unsupported number of gpus and mpisize. Should be 1 or a multiple of n_gpus=%d\n", n_gpus);
        }

        return d;
    }

    void init_gpu(device & d)
    {
        CHECK(hipSetDevice(d->gpu_idx));
        CHECK(hipSetDeviceFlags(hipDeviceScheduleYield));
        CHECK(hipFree(nullptr));
        CHECK(hipDeviceSynchronize());
    }

    void set_device(device & d)
    {
        #pragma omp parallel
        {
            CHECK(hipSetDevice(d->gpu_idx));
        }
    }

    void queue_create(queue & q)
    {
        q = std::make_shared<_queue>();
        CHECK(hipStreamCreate(&q->stream));
    }
    void queue_destroy(queue & q)
    {
        CHECK(hipStreamDestroy(q->stream));
        q.reset();
    }

    void queue_async_barrier(const std::vector<queue> & waitfor, const std::vector<queue> & waitin)
    {
        // made mainly for 1:N or N:1 scenarios
        std::vector<hipEvent_t> events(waitfor.size());
        for(size_t i = 0; i < waitfor.size(); i++) {
            CHECK(hipEventCreate(&events[i]));
            CHECK(hipEventRecord(events[i], waitfor[i]->stream));
        }
        for(size_t j = 0; j < waitin.size(); j++) {
            for(size_t k = 0; k < events.size(); k++) {
                CHECK(hipStreamWaitEvent(waitin[j]->stream, events[k], 0));
            }
        }
        for(size_t k = 0; k < events.size(); k++) {
            CHECK(hipEventDestroy(events[k]));
        }
    }

    void queue_wait(queue & q)
    {
        CHECK(hipStreamSynchronize(q->stream));
    }

    void device_wait()
    {
        CHECK(hipDeviceSynchronize());
    }

    size_t get_device_memory_capacity()
    {
        hipDeviceProp_t props;
        int gpu_idx;
        CHECK(hipGetDevice(&gpu_idx));
        CHECK(hipGetDeviceProperties(&props, gpu_idx));
        return props.totalGlobalMem;
    }

    size_t get_device_memory_free()
    {
        size_t size_free;
        size_t size_total;
        CHECK(hipMemGetInfo(&size_free, &size_total));
        return size_free;
    }

    size_t get_natural_pitch_align()
    {
        return 512;
    }

    void * memalloc_device(size_t num_bytes)
    {
        void * ptr;
        CHECK(hipMalloc(&ptr, num_bytes));
        return ptr;
    }

    void * memalloc_device_2d(size_t num_chunks, size_t bytes_per_chunk, size_t & pitch)
    {
        void * ptr;
        CHECK(hipMallocPitch(&ptr, &pitch, bytes_per_chunk, num_chunks));
        return ptr;
    }

    void memfree_device(void * ptr)
    {
        CHECK(hipFree(ptr));
    }

    void memalloc_device_max(void * & memory, size_t & memory_size_B, size_t max_needed)
    {
        size_t keep_free_percent = 5;
        size_t size_free;
        size_t size_total;
        CHECK(hipMemGetInfo(&size_free, &size_total));
        size_t can_allocate_max = ((100 - keep_free_percent) * size_free) / 100;
        memory_size_B = std::min(can_allocate_max, max_needed);
        CHECK(hipMalloc(&memory, memory_size_B));
    }

    void * memalloc_hostpinned(size_t num_bytes)
    {
        void * ptr;
        CHECK(hipHostMalloc(&ptr, num_bytes));
        return ptr;
    }

    void memfree_hostpinned(void * ptr)
    {
        CHECK(hipHostFree(ptr));
    }

    void submit_host_function(queue & q, const std::function<void(void)> & f)
    {
        std::function<void(void)> * func_ptr = new std::function<void(void)>(f);

        CHECK(hipStreamAddCallback(q->stream, [](hipStream_t /*str*/, hipError_t /*err*/, void * arg){
            std::function<void(void)> * func_ptr = reinterpret_cast<std::function<void(void)>*>(arg);
            (*func_ptr)();
            delete func_ptr;
        }, func_ptr, 0));

        // hipLaunchHostFunc is still marked as beta in rocm-6.0.0
    }

    template<typename T>
    void copy_submit(queue & q, T * dst, T const * src, size_t num_elements)
    {
        size_t sizeof_T;
        if constexpr(std::is_same_v<T,void>) sizeof_T = 1;
        else sizeof_T = sizeof(T);
        size_t size = num_elements * sizeof_T;
        CHECK(hipMemcpyAsync(dst, src, size, hipMemcpyDefault, q->stream));
    }

    template<typename T, typename I, typename Ao, typename Ai>
    void copy_submit(queue & q, Vector_Dense<T,I,Ao> & output, const Vector_Dense<T,I,Ai> & input)
    {
        if(output.size != input.size) eslog::error("copy submit: output vector has wrong dimensions\n");
        copy_submit(q, output.vals, input.vals, input.size);
    }

    template<typename T, typename I, typename Ao, typename Ai>
    void copy_submit(queue & q, Matrix_Dense<T,I,Ao> & output, const Matrix_Dense<T,I,Ai> & input)
    {
        if(output.nrows != input.nrows || output.ncols != input.ncols) eslog::error("copy submit: output matrix has wrong dimensions\n");
        CHECK(hipMemcpy2DAsync(output.vals, output.get_ld() * sizeof(T), input.vals, input.get_ld() * sizeof(T), input.ncols * sizeof(T), input.nrows, hipMemcpyDefault, q->stream));
    }

    template<typename T, typename I, typename Ao, typename Ai>
    void copy_submit(queue & q, Matrix_CSR<T,I,Ao> & output, const Matrix_CSR<T,I,Ai> & input, bool copy_pattern, bool copy_vals)
    {
        if(output.nrows != input.nrows || output.ncols != input.ncols || output.nnz != input.nnz) eslog::error("copy submit: output matrix has wrong dimensions\n");
        if(copy_pattern) copy_submit(q, output.rows, input.rows, input.nrows+1);
        if(copy_pattern) copy_submit(q, output.cols, input.cols, input.nnz);
        if(copy_vals)    copy_submit(q, output.vals, input.vals, input.nnz);
    }

    template<typename T>
    void copy_submit(queue & q, PermutationView_new<T> & src, PermutationView_new<T> & dst)
    {
        if(src.size != dst.size) eslog::error("copy submit: output permutation has wrong dimensions\n");
        copy_submit(q, dst.src_to_dst, dst.src_to_dst, input.size);
        copy_submit(q, dst.dst_to_src, dst.dst_to_src, input.size);
    }

    template<typename T>
    void copy_submit(queue & /*q*/, VectorDenseView_new<T> & /*src*/, VectorDenseView_new<T> & /*dst*/)
    {
        eslog::error("not supported yet\n");
    }

    template<typename T>
    void copy_submit(queue & /*q*/, MatrixDenseView_new<T> & /*src*/, MatrixDenseView_new<T> & /*dst*/)
    {
        eslog::error("not supported yet\n");
    }

    template<typename T, typename I>
    void copy_submit(queue & /*q*/, MatrixCsxView_new<T,I> & /*src*/, MatrixCsxView_new<T,I> & /*dst*/, bool /*copy_pattern*/, bool /*copy_vals*/)
    {
        eslog::error("not supported yet\n");
    }

    void memset_submit(queue & q, void * ptr, size_t num_bytes, char val)
    {
        CHECK(hipMemsetAsync(ptr, val, num_bytes, q->stream));
    }

}
}
}

#include "gpu/gpu_management.inst.hpp"

#endif
#endif
