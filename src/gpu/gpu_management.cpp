
#if !defined(HAVE_CUDA) && !defined(HAVE_ROCM)

#include "gpu_management.h"

namespace espreso {
namespace gpu {
namespace mgm {

    struct _device {};

    struct _queue {};

    void* Ad::allocate(size_t num_bytes)
    {
        return nullptr;
    }

    void Ad::deallocate(void *ptr)
    {

    }

    void* Ah::allocate(size_t num_bytes)
    {
        return nullptr;
    }

    void Ah::deallocate(void *ptr)
    {

    }

    device get_device_by_mpi(int mpi_rank, int mpi_size) { return device{}; }

    void init_gpu(device & d) {}

    void set_device(const device & d) {}

    void queue_create(queue & q, device & d) {}

    void queue_destroy(queue & q) {}

    void queue_async_barrier(const std::vector<queue> & waitfor, const std::vector<queue> & waitin) {}

    void queue_wait(queue & q) {}

    void device_wait(device & d) {}

    size_t get_device_memory_capacity(device & d) { return 0; }

    void * memalloc_device(device & d, size_t num_bytes) { return nullptr; }

    void memfree_device(device & d, void * ptr) {}

    void memalloc_device_max(device & d, void * & memory, size_t & memory_size_B, size_t max_needed) {}

    void submit_host_function(queue & q, const std::function<void(void)> & c) {}

    template<typename T, typename I>
    void copy_submit_h2d(queue & q, T * dst, const T * src, I num_elements) {}

    template<typename T, typename I>
    void copy_submit_d2h(queue & q, T * dst, const T * src, I num_elements) {}

    template<typename T, typename I, typename Ao, typename Ai>
    void copy_submit_d2h(queue & q, Vector_Dense<T,I,Ao> & output, const Vector_Dense<T,I,Ai> & input) {}

    template<typename T, typename I, typename Ao, typename Ai>
    void copy_submit_h2d(queue & q, Vector_Dense<T,I,Ao> & output, const Vector_Dense<T,I,Ai> & input) {}

    template<typename T, typename I, typename Ao, typename Ai>
    void copy_submit_d2h(queue & q, Matrix_Dense<T,I,Ao> & output, const Matrix_Dense<T,I,Ai> & input) {}

    template<typename T, typename I, typename Ao, typename Ai>
    void copy_submit_h2d(queue & q, Matrix_Dense<T,I,Ao> & output, const Matrix_Dense<T,I,Ai> & input) {}
    
    template<typename T, typename I, typename Ao, typename Ai>
    void copy_submit_d2h(queue & q, Matrix_CSR<T,I,Ao> & output, const Matrix_CSR<T,I,Ai> & input, bool copy_pattern = true, bool copy_vals = true) {}

    template<typename T, typename I, typename Ao, typename Ai>
    void copy_submit_h2d(queue & q, Matrix_CSR<T,I,Ao> & output, const Matrix_CSR<T,I,Ai> & input, bool copy_pattern = true, bool copy_vals = true) {}

    void memset_submit(queue & q, void * ptr, size_t num_bytes, char val) {}



    #define INSTANTIATE(T,I,Ahost,Adevice) \
    template void copy_submit_h2d<T,I,Adevice,Ahost  >(queue & q, T * dst, const T * src, I num_elements); \
    template void copy_submit_d2h<T,I,Ahost,  Adevice>(queue & q, T * dst, const T * src, I num_elements); \
    template void copy_submit_h2d<T,I,Adevice,Ahost  >(queue & q, Vector_Dense<T,I,Ao> & output, const Vector_Dense<T,I,Ai> & input); \
    template void copy_submit_d2h<T,I,Ahost,  Adevice>(queue & q, Vector_Dense<T,I,Ao> & output, const Vector_Dense<T,I,Ai> & input); \
    template void copy_submit_h2d<T,I,Adevice,Ahost  >(queue & q, Matrix_Dense<T,I,Ao> & output, const Matrix_Dense<T,I,Ai> & input); \
    template void copy_submit_d2h<T,I,Ahost,  Adevice>(queue & q, Matrix_Dense<T,I,Ao> & output, const Matrix_Dense<T,I,Ai> & input); \
    template void copy_submit_h2d<T,I,Adevice,Ahost  >(queue & q, Matrix_CSR<T,I,Ao> & output, const Matrix_CSR<T,I,Ai> & input, bool copy_pattern = true, bool copy_vals = true); \
    template void copy_submit_d2h<T,I,Ahost,  Adevice>(queue & q, Matrix_CSR<T,I,Ao> & output, const Matrix_CSR<T,I,Ai> & input, bool copy_pattern = true, bool copy_vals = true);
        // INSTANTIATE(float,                int32_t, mgm::Ah,       mgm::Ad)
        INSTANTIATE(double,               int32_t, mgm::Ah,       mgm::Ad)
        // INSTANTIATE(std::complex<float >, int32_t, mgm::Ah,       mgm::Ad)
        // INSTANTIATE(std::complex<double>, int32_t, mgm::Ah,       mgm::Ad)
        // INSTANTIATE(float,                int64_t, mgm::Ah,       mgm::Ad)
        // INSTANTIATE(double,               int64_t, mgm::Ah,       mgm::Ad)
        // INSTANTIATE(std::complex<float >, int64_t, mgm::Ah,       mgm::Ad)
        // INSTANTIATE(std::complex<double>, int64_t, mgm::Ah,       mgm::Ad)
        // INSTANTIATE(float,                int32_t, mgm::Ah,       cbmba_d)
        // INSTANTIATE(double,               int32_t, mgm::Ah,       cbmba_d)
        // INSTANTIATE(std::complex<float >, int32_t, mgm::Ah,       cbmba_d)
        // INSTANTIATE(std::complex<double>, int32_t, mgm::Ah,       cbmba_d)
        // INSTANTIATE(float,                int64_t, mgm::Ah,       cbmba_d)
        // INSTANTIATE(double,               int64_t, mgm::Ah,       cbmba_d)
        // INSTANTIATE(std::complex<float >, int64_t, mgm::Ah,       cbmba_d)
        // INSTANTIATE(std::complex<double>, int64_t, mgm::Ah,       cbmba_d)
        // INSTANTIATE(float,                int32_t, cpu_allocator, mgm::Ad)
        // INSTANTIATE(double,               int32_t, cpu_allocator, mgm::Ad)
        // INSTANTIATE(std::complex<float >, int32_t, cpu_allocator, mgm::Ad)
        // INSTANTIATE(std::complex<double>, int32_t, cpu_allocator, mgm::Ad)
        // INSTANTIATE(float,                int64_t, cpu_allocator, mgm::Ad)
        // INSTANTIATE(double,               int64_t, cpu_allocator, mgm::Ad)
        // INSTANTIATE(std::complex<float >, int64_t, cpu_allocator, mgm::Ad)
        // INSTANTIATE(std::complex<double>, int64_t, cpu_allocator, mgm::Ad)
        // INSTANTIATE(float,                int32_t, cpu_allocator, cbmba_d)
        // INSTANTIATE(double,               int32_t, cpu_allocator, cbmba_d)
        // INSTANTIATE(std::complex<float >, int32_t, cpu_allocator, cbmba_d)
        // INSTANTIATE(std::complex<double>, int32_t, cpu_allocator, cbmba_d)
        // INSTANTIATE(float,                int64_t, cpu_allocator, cbmba_d)
        // INSTANTIATE(double,               int64_t, cpu_allocator, cbmba_d)
        // INSTANTIATE(std::complex<float >, int64_t, cpu_allocator, cbmba_d)
        // INSTANTIATE(std::complex<double>, int64_t, cpu_allocator, cbmba_d)
    #undef INSTANTIATE

}
}
}

#endif
