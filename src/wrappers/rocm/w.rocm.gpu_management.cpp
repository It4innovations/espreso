
#ifdef HAVE_ROCM

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

    namespace
    {
        template<typename T, typename I>
        static void _copy_submit(T * dst, const T * src, I num_elements, hipMemcpyKind dir, hipStream_t stream)
        {
            size_t sizeof_T;
            if constexpr(std::is_same_v<T,void>) sizeof_T = 1;
            else sizeof_T = sizeof(T);
            size_t size = num_elements * sizeof_T;
            CHECK(hipMemcpyAsync(dst, src, size, dir, stream));
        }

        template<typename T, typename I, typename A1, typename A2>
        static void _copy_submit(Vector_Dense<T,I,A1> & output, const Vector_Dense<T,I,A2> & input, hipMemcpyKind direction, hipStream_t stream)
        {
            if(output.size != input.size) eslog::error("copy submit: output vector has wrong dimensions\n");
            CHECK(hipMemcpyAsync(output.vals, input.vals, input.size * sizeof(T), direction, stream));
        }

        template<typename T, typename I, typename A1, typename A2>
        static void _copy_submit(Matrix_Dense<T,I,A1> & output, const Matrix_Dense<T,I,A2> & input, hipMemcpyKind direction, hipStream_t stream)
        {
            if(output.nrows != input.nrows || output.ncols != input.ncols) eslog::error("copy submit: output matrix has wrong dimensions\n");
            CHECK(hipMemcpy2DAsync(output.vals, output.get_ld() * sizeof(T), input.vals, input.get_ld() * sizeof(T), input.ncols * sizeof(T), input.nrows, direction, stream));
        }

        template<typename T, typename I, typename A1, typename A2>
        static void _copy_submit(Matrix_CSR<T,I,A1> & output, const Matrix_CSR<T,I,A2> & input, hipMemcpyKind direction, hipStream_t stream, bool copy_pattern, bool copy_vals)
        {
            if(output.nrows != input.nrows || output.ncols != input.ncols || output.nnz != input.nnz) eslog::error("copy submit: output matrix has wrong dimensions\n");
            if(copy_pattern) CHECK(hipMemcpyAsync(output.rows, input.rows, (input.nrows+1) * sizeof(I), direction, stream));
            if(copy_pattern) CHECK(hipMemcpyAsync(output.cols, input.cols,  input.nnz      * sizeof(I), direction, stream));
            if(copy_vals)    CHECK(hipMemcpyAsync(output.vals, input.vals,  input.nnz      * sizeof(T), direction, stream));
        }
    }

    device get_device_by_mpi(int mpi_rank, int mpi_size)
    {
        // static constexpr int rank_gpu_map[] = {4,5,2,3,6,7,0,1}; // assuming LUMI-G
        // on LUMI, when using e.g. `salloc --ntasks=8 --gpus-per-task=1`, each task (=rank) can see only a single gpu
        device d = std::make_shared<_device>();
        d->gpu_idx = 0;
        return d;
    }

    void init_gpu(device & d)
    {
        CHECK(hipSetDevice(d->gpu_idx));
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
        for(size_t i = 0; i < waitfor.size(); i++)
        {
            CHECK(hipEventCreate(&events[i]));
            CHECK(hipEventRecord(events[i], waitfor[i]->stream));
        }
        for(size_t j = 0; j < waitin.size(); j++)
        {
            for(size_t k = 0; k < events.size(); k++)
            {
                CHECK(hipStreamWaitEvent(waitin[j]->stream, events[k], 0));
            }
        }
        for(size_t k = 0; k < events.size(); k++)
        {
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

    void * memalloc_device(size_t num_bytes)
    {
        void * ptr;
        CHECK(hipMalloc(&ptr, num_bytes));
        return ptr;
    }

    void memfree_device(void * ptr)
    {
        CHECK(hipFree(ptr));
    }

    void memalloc_device_max(void * & memory, size_t & memory_size_B, size_t max_needed)
    {
        size_t coef_percent = 95;

        memory_size_B = std::min(max_needed, get_device_memory_capacity());
        while(memory_size_B > 0)
        {
            hipError_t err = hipMalloc(&memory, memory_size_B);

            if(err == hipSuccess)
            {
                err = hipGetLastError(); // reset the internal last error variable
                return;
            }
            if(err != hipErrorMemoryAllocation) CHECK(err);

            memory_size_B = (memory_size_B * coef_percent) / 100;
        }

        eslog::error("could not allocate any gpu memory\n");
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

    template<typename T, typename I>
    void copy_submit_h2d(queue & q, T * dst, T const * src, I num_elements)
    {
        _copy_submit(dst, src, num_elements, hipMemcpyHostToDevice, q->stream);
    }
    template<typename T, typename I>
    void copy_submit_d2h(queue & q, T * dst, T const * src, I num_elements)
    {
        _copy_submit(dst, src, num_elements, hipMemcpyDeviceToHost, q->stream);
    }

    template<typename T, typename I, typename Ao, typename Ai>
    void copy_submit_d2h(queue & q, Vector_Dense<T,I,Ao> & output, const Vector_Dense<T,I,Ai> & input)
    {
        static_assert(Ao::is_data_host_accessible, "output vector data has to be host accessible");
        static_assert(Ai::is_data_device_accessible, "input vector data has to be device accessible");
        _copy_submit(output, input, hipMemcpyDeviceToHost, q->stream);
    }

    template<typename T, typename I, typename Ao, typename Ai>
    void copy_submit_h2d(queue & q, Vector_Dense<T,I,Ao> & output, const Vector_Dense<T,I,Ai> & input)
    {
        static_assert(Ao::is_data_device_accessible, "output vector data has to be device accessible");
        static_assert(Ai::is_data_host_accessible, "input vector data has to be host accessible");
        _copy_submit(output, input, hipMemcpyHostToDevice, q->stream);
    }

    template<typename T, typename I, typename Ao, typename Ai>
    void copy_submit_d2h(queue & q, Matrix_Dense<T,I,Ao> & output, const Matrix_Dense<T,I,Ai> & input)
    {
        static_assert(Ao::is_data_host_accessible, "output matrix data has to be host accessible");
        static_assert(Ai::is_data_device_accessible, "input matrix data has to be device accessible");
        _copy_submit(output, input, hipMemcpyDeviceToHost, q->stream);
    }

    template<typename T, typename I, typename Ao, typename Ai>
    void copy_submit_h2d(queue & q, Matrix_Dense<T,I,Ao> & output, const Matrix_Dense<T,I,Ai> & input)
    {
        static_assert(Ao::is_data_device_accessible, "output matrix data has to be device accessible");
        static_assert(Ai::is_data_host_accessible, "input matrix data has to be host accessible");
        _copy_submit(output, input, hipMemcpyHostToDevice, q->stream);
    }

    template<typename T, typename I, typename Ao, typename Ai>
    void copy_submit_d2h(queue & q, Matrix_CSR<T,I,Ao> & output, const Matrix_CSR<T,I,Ai> & input, bool copy_pattern, bool copy_vals)
    {
        static_assert(Ao::is_data_host_accessible, "output matrix data has to be host accessible");
        static_assert(Ai::is_data_device_accessible, "input matrix data has to be device accessible");
        _copy_submit(output, input, hipMemcpyDeviceToHost, q->stream, copy_pattern, copy_vals);
    }

    template<typename T, typename I, typename Ao, typename Ai>
    void copy_submit_h2d(queue & q, Matrix_CSR<T,I,Ao> & output, const Matrix_CSR<T,I,Ai> & input, bool copy_pattern, bool copy_vals)
    {
        static_assert(Ao::is_data_device_accessible, "output matrix data has to be device accessible");
        static_assert(Ai::is_data_host_accessible, "input matrix data has to be host accessible");
        _copy_submit(output, input, hipMemcpyHostToDevice, q->stream, copy_pattern, copy_vals);
    }

    void memset_submit(queue & q, void * ptr, size_t num_bytes, char val)
    {
        CHECK(hipMemsetAsync(ptr, val, num_bytes, q->stream));
    }

}
}
}

#include "gpu/gpu_management_inst.hpp"

#endif
