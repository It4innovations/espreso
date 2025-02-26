
#ifdef HAVE_CUDA
#ifdef ESPRESO_USE_WRAPPER_GPU_CUDA

#include "gpu/gpu_management.h"
#include "w.cuda.gpu_management.h"
#include "basis/utilities/cbmb_allocator.h"

#include <omp.h>
#include <complex>



namespace espreso {
namespace gpu {
namespace mgm {

    gpu_wrapper_impl get_implementation()
    {
        return gpu_wrapper_impl::CUDA;
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
            // assuming that there are as many ranks on a node as gpus
            int local_node_rank = mpi_rank % n_gpus;
            d->gpu_idx = rank_gpu_map[local_node_rank];
        }
        else {
            eslog::error("unsupported number of gpus and mpisize. Should be 1 or a multiple of n_gpus=%d\n", n_gpus);
        }

        return d;
    }

    void init_gpu(device & d)
    {
        CHECK(cudaSetDevice(d->gpu_idx));
        CHECK(cudaSetDeviceFlags(cudaDeviceScheduleYield));
        CHECK(cudaFree(nullptr));
        CHECK(cudaDeviceSynchronize());
    }

    void set_device(device & d)
    {
        #pragma omp parallel
        {
            CHECK(cudaSetDevice(d->gpu_idx));
        }
    }

    void queue_create(queue & q)
    {
        q = std::make_shared<_queue>();
        CHECK(cudaStreamCreate(&q->stream));
    }

    void queue_destroy(queue & q)
    {
        CHECK(cudaStreamDestroy(q->stream));
        q.reset();
    }

    void queue_async_barrier(const std::vector<queue> & waitfor, const std::vector<queue> & waitin)
    {
        // made mainly for 1:N or N:1 scenarios
        std::vector<cudaEvent_t> events(waitfor.size());
        for(size_t i = 0; i < waitfor.size(); i++) {
            CHECK(cudaEventCreate(&events[i]));
            CHECK(cudaEventRecord(events[i], waitfor[i]->stream));
        }
        for(size_t j = 0; j < waitin.size(); j++) {
            for(size_t k = 0; k < events.size(); k++) {
                CHECK(cudaStreamWaitEvent(waitin[j]->stream, events[k], 0));
            }
        }
        for(size_t k = 0; k < events.size(); k++) {
            CHECK(cudaEventDestroy(events[k]));
        }
    }

    void queue_wait(queue & q)
    {
        CHECK(cudaStreamSynchronize(q->stream));
    }

    void device_wait()
    {
        CHECK(cudaDeviceSynchronize());
    }

    size_t get_device_memory_capacity()
    {
        cudaDeviceProp props;
        int gpu_idx;
        CHECK(cudaGetDevice(&gpu_idx));
        CHECK(cudaGetDeviceProperties(&props, gpu_idx));
        return props.totalGlobalMem;
    }

    size_t get_device_memory_free()
    {
        size_t size_free;
        size_t size_total;
        CHECK(cudaMemGetInfo(&size_free, &size_total));
        return size_free;
    }

    void * memalloc_device(size_t num_bytes)
    {
        void * ptr;
        CHECK(cudaMalloc(&ptr, num_bytes));
        return ptr;
    }

    void memfree_device(void * ptr)
    {
        CHECK(cudaFree(ptr));
    }

    void memalloc_device_max(void * & memory, size_t & memory_size_B, size_t max_needed)
    {
        size_t keep_free_percent = 5;
        size_t size_free;
        size_t size_total;
        CHECK(cudaMemGetInfo(&size_free, &size_total));
        size_t can_allocate_max = ((100 - keep_free_percent) * size_free) / 100;
        memory_size_B = std::min(can_allocate_max, max_needed);
        CHECK(cudaMalloc(&memory, memory_size_B));
    }

    void * memalloc_hostpinned(size_t num_bytes)
    {
        void * ptr;
        CHECK(cudaMallocHost(&ptr, num_bytes));
        return ptr;
    }

    void memfree_hostpinned(void * ptr)
    {
        CHECK(cudaFreeHost(ptr));
    }

    void submit_host_function(queue & q, const std::function<void(void)> & f)
    {
        std::function<void(void)> * func_ptr = new std::function<void(void)>(f);

        CHECK(cudaLaunchHostFunc(q->stream, [](void * arg){
            std::function<void(void)> * func_ptr = reinterpret_cast<std::function<void(void)>*>(arg);
            (*func_ptr)();
            delete func_ptr;
        }, func_ptr));
    }

    template<typename T>
    void copy_submit(queue & q, T * dst, T const * src, size_t num_elements)
    {
        size_t sizeof_T;
        if constexpr(std::is_same_v<T,void>) sizeof_T = 1;
        else sizeof_T = sizeof(T);
        size_t size = num_elements * sizeof_T;
        CHECK(cudaMemcpyAsync(dst, src, size, cudaMemcpyDefault, q->stream));
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
        CHECK(cudaMemcpy2DAsync(output.vals, output.get_ld() * sizeof(T), input.vals, input.get_ld() * sizeof(T), input.ncols * sizeof(T), input.nrows, cudaMemcpyDefault, q->stream));
    }

    template<typename T, typename I, typename Ao, typename Ai>
    void copy_submit(queue & q, Matrix_CSR<T,I,Ao> & output, const Matrix_CSR<T,I,Ai> & input, bool copy_pattern, bool copy_vals)
    {
        if(output.nrows != input.nrows || output.ncols != input.ncols || output.nnz != input.nnz) eslog::error("copy submit: output matrix has wrong dimensions\n");
        if(copy_pattern) copy_submit(q, output.rows, input.rows, input.nrows+1);
        if(copy_pattern) copy_submit(q, output.cols, input.cols, input.nnz);
        if(copy_vals)    copy_submit(q, output.vals, input.vals, input.nnz);
    }

    void memset_submit(queue & q, void * ptr, size_t num_bytes, char val)
    {
        CHECK(cudaMemsetAsync(ptr, val, num_bytes, q->stream));
    }

}
}
}

#include "gpu/gpu_management.inst.hpp"

#endif
#endif
