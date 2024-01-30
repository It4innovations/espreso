
#ifdef HAVE_CUDA

#include "gpu/gpu_management.h"

#include <omp.h>
#include "w.cuda.common.h"



namespace espreso {
namespace gpu {
namespace mgm {

    namespace
    {
        template<typename T, typename I>
        static void _copy_submit(T * dst, const T * src, I num_elements, cudaMemcpyKind dir, cudaStream_t stream)
        {
            size_t sizeof_T;
            if constexpr(std::is_same_v<T,void>) sizeof_T = 1;
            else sizeof_T = sizeof(T);
            size_t size = num_elements * sizeof_T;
            CHECK(cudaMemcpyAsync(dst, src, size, dir, stream));
        }

        template<typename T, typename I, typename A1, typename A2>
        static void _copy_submit(Vector_Dense<T,I,A1> & output, const Vector_Dense<T,I,A2> & input, cudaMemcpyKind direction, cudaStream_t stream)
        {
            if(output.size != input.size) eslog::error("copy submit: output vector has wrong dimensions");
            CHECK(cudaMemcpyAsync(output.vals, input.vals, input.size * sizeof(T), direction, stream));
        }

        template<typename T, typename I, typename A1, typename A2>
        static void _copy_submit(Matrix_Dense<T,I,A1> & output, const Matrix_Dense<T,I,A2> & input, cudaMemcpyKind direction, cudaStream_t stream)
        {
            if(output.nrows != input.nrows || output.ncols != input.ncols) eslog::error("copy submit: output matrix has wrong dimensions");
            CHECK(cudaMemcpy2DAsync(output.vals, output.get_ld() * sizeof(T), input.vals, input.get_ld() * sizeof(T), input.ncols * sizeof(T), input.nrows, direction, stream));
        }

        template<typename T, typename I, typename A1, typename A2>
        static void _copy_submit(Matrix_CSR<T,I,A1> & output, const Matrix_CSR<T,I,A2> & input, cudaMemcpyKind direction, cudaStream_t stream, bool copy_pattern, bool copy_vals)
        {
            if(output.nrows != input.nrows || output.ncols != input.ncols || output.nnz != input.nnz) eslog::error("copy submit: output matrix has wrong dimensions");
            if(copy_pattern) CHECK(cudaMemcpyAsync(output.rows, input.rows, (input.nrows+1) * sizeof(I), direction, stream));
            if(copy_pattern) CHECK(cudaMemcpyAsync(output.cols, input.cols,  input.nnz      * sizeof(I), direction, stream));
            if(copy_vals)    CHECK(cudaMemcpyAsync(output.vals, input.vals,  input.nnz      * sizeof(T), direction, stream));
        }
    }

    struct device
    {
        int gpu_idx;
    };

    struct queue
    {
        cudaStream_t stream;
    };

    inline void * Ad::allocate(size_t num_bytes)
    {
        void * ptr;
        CHECK(cudaMalloc(&ptr, num_bytes));
        return ptr;
    }

    template<typename T>
    inline T * Ad::allocate(size_t count)
    {
        return reinterpret_cast<T*>(allocate(count * sizeof(T)));
    }

    template<typename T>
    inline void Ad::deallocate(T * ptr)
    { 
        CHECK(cudaFree(ptr));
    }

    inline void * Ah::allocate(size_t num_bytes)
    {
        void * ptr;
        CHECK(cudaMallocHost(&ptr, num_bytes));
        return ptr;
    }

    template<typename T>
    inline T * Ah::allocate(size_t count)
    {
        return reinterpret_cast<T*>(allocate(count * sizeof(T)));
    }

    template<typename T>
    inline void Ah::deallocate(T * ptr)
    { 
        CHECK(cudaFreeHost(ptr));
    }

    static device get_device_by_mpi(int mpi_rank, int mpi_size)
    {
        static constexpr int rank_gpu_map[] = {2,3,0,1,6,7,4,5}; // assuming Karolina GPU
        device d;
        if(mpi_size == 1) {
            d.gpu_idx = 0;
        }
        else if(mpi_size % 8 == 0) {
            int local_node_rank = mpi_rank % 8;
            d.gpu_idx = rank_gpu_map[local_node_rank];
        }
        else {
            eslog::error("unsupported number of gpus and mpisize. Should be 1 or a multiple of 8\n");
            d.gpu_idx = -1;
        }

        return d;
    }

    static void init_gpu(const device & d)
    {
        CHECK(cudaSetDevice(d.gpu_idx));
        CHECK(cudaFree(nullptr));
        CHECK(cudaDeviceSynchronize());
    }

    static void set_device(const device & d)
    {
        #pragma omp parallel
        {
            CHECK(cudaSetDevice(d.gpu_idx));
        }
    }

    static void queue_create(queue & q, device & /*d*/)
    {
        CHECK(cudaStreamCreate(&q.stream));
    }

    static void queue_destroy(queue & q)
    {
        CHECK(cudaStreamDestroy(q.stream));
    }

    static void queue_async_barrier(const std::vector<queue> & waitfor, const std::vector<queue> & waitin)
    {
        // made mainly for 1:N or N:1 scenarios
        std::vector<cudaEvent_t> events(waitfor.size());
        for(size_t i = 0; i < waitfor.size(); i++)
        {
            CHECK(cudaEventCreate(&events[i]));
            CHECK(cudaEventRecord(events[i], waitfor[i].stream));
        }
        for(size_t j = 0; j < waitin.size(); j++)
        {
            for(size_t k = 0; k < events.size(); k++)
            {
                CHECK(cudaStreamWaitEvent(waitin[j].stream, events[k], 0));
            }
        }
        for(size_t k = 0; k < events.size(); k++)
        {
            CHECK(cudaEventDestroy(events[k]));
        }
    }

    static void queue_wait(queue & q)
    {
        CHECK(cudaStreamSynchronize(q.stream));
    }

    static void device_wait(device & /*d*/)
    {
        CHECK(cudaDeviceSynchronize());
    }

    static size_t get_device_memory_capacity(device & d)
    {
        cudaDeviceProp props;
        CHECK(cudaGetDeviceProperties(&props, d.gpu_idx));
        return props.totalGlobalMem;
    }

    static void * memalloc_device(device & /*d*/, size_t num_bytes)
    {
        void * ptr;
        CHECK(cudaMalloc(&ptr, num_bytes));
        return ptr;
    }

    template<typename T>
    static void memfree_device(device & /*d*/, T * ptr)
    {
        CHECK(cudaFree(ptr));
    }

    static void memalloc_device_max(device & d, void * & memory, size_t & memory_size_B, size_t max_needed)
    {
        size_t coef_percent = 95;

        memory_size_B = std::min(max_needed, get_device_memory_capacity(d));
        while(memory_size_B > 0)
        {
            cudaMalloc(&memory, memory_size_B);
            cudaError_t err = cudaGetLastError();

            if(err == cudaSuccess) return;
            if(err != cudaErrorMemoryAllocation) CHECK(err);

            memory_size_B = (memory_size_B * coef_percent) / 100;
        }

        eslog::error("could not allocate any gpu memory");
    }


    template<typename C>
    static void submit_host_function(queue & q, C && c)
    {
        C * cp = new C(std::move(c));

        CHECK(cudaLaunchHostFunc(q.stream, [](void * arg){
            C * cpi = reinterpret_cast<C*>(arg);
            (*cpi)();
            delete cpi;
        }, cp));
    }

    template<typename T, typename I>
    static void copy_submit_h2d(queue & q, T * dst, const T * src, I num_elements)
    {
        _copy_submit(dst, src, num_elements, cudaMemcpyHostToDevice, q.stream);
    }

    template<typename T, typename I>
    static void copy_submit_d2h(queue & q, T * dst, const T * src, I num_elements)
    {
        _copy_submit(dst, src, num_elements, cudaMemcpyDeviceToHost, q.stream);
    }

    template<typename T, typename I, typename Ao, typename Ai>
    static void copy_submit_d2h(queue & q, Vector_Dense<T,I,Ao> & output, const Vector_Dense<T,I,Ai> & input)
    {
        static_assert(Ao::is_data_host_accessible, "output vector data has to be host accessible");
        static_assert(Ai::is_data_device_accessible, "input vector data has to be device accessible");
        _copy_submit(output, input, cudaMemcpyDeviceToHost, q.stream);
    }

    template<typename T, typename I, typename Ao, typename Ai>
    static void copy_submit_h2d(queue & q, Vector_Dense<T,I,Ao> & output, const Vector_Dense<T,I,Ai> & input)
    {
        static_assert(Ao::is_data_device_accessible, "output vector data has to be device accessible");
        static_assert(Ai::is_data_host_accessible, "input vector data has to be host accessible");
        _copy_submit(output, input, cudaMemcpyHostToDevice, q.stream);
    }

    template<typename T, typename I, typename Ao, typename Ai>
    static void copy_submit_d2h(queue & q, Matrix_Dense<T,I,Ao> & output, const Matrix_Dense<T,I,Ai> & input)
    {
        static_assert(Ao::is_data_host_accessible, "output matrix data has to be host accessible");
        static_assert(Ai::is_data_device_accessible, "input matrix data has to be device accessible");
        _copy_submit(output, input, cudaMemcpyDeviceToHost, q.stream);
    }

    template<typename T, typename I, typename Ao, typename Ai>
    static void copy_submit_h2d(queue & q, Matrix_Dense<T,I,Ao> & output, const Matrix_Dense<T,I,Ai> & input)
    {
        static_assert(Ao::is_data_device_accessible, "output matrix data has to be device accessible");
        static_assert(Ai::is_data_host_accessible, "input matrix data has to be host accessible");
        _copy_submit(output, input, cudaMemcpyHostToDevice, q.stream);
    }

    template<typename T, typename I, typename Ao, typename Ai>
    static void copy_submit_d2h(queue & q, Matrix_CSR<T,I,Ao> & output, const Matrix_CSR<T,I,Ai> & input, bool copy_pattern, bool copy_vals)
    {
        static_assert(Ao::is_data_host_accessible, "output matrix data has to be host accessible");
        static_assert(Ai::is_data_device_accessible, "input matrix data has to be device accessible");
        _copy_submit(output, input, cudaMemcpyDeviceToHost, q.stream, copy_pattern, copy_vals);
    }

    template<typename T, typename I, typename Ao, typename Ai>
    static void copy_submit_h2d(queue & q, Matrix_CSR<T,I,Ao> & output, const Matrix_CSR<T,I,Ai> & input, bool copy_pattern, bool copy_vals)
    {
        static_assert(Ao::is_data_device_accessible, "output matrix data has to be device accessible");
        static_assert(Ai::is_data_host_accessible, "input matrix data has to be host accessible");
        _copy_submit(output, input, cudaMemcpyHostToDevice, q.stream, copy_pattern, copy_vals);
    }

    template<typename T, typename I>
    static void memset_submit(queue & q, T * ptr, I num_bytes, char val)
    {
        CHECK(cudaMemsetAsync(ptr, val, num_bytes, q.stream));
    }

}
}
}

#endif
