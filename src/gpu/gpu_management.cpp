
#if !defined(HAVE_CUDA) && !defined(HAVE_ROCM)

#include "gpu_management.h"
#include "basis/utilities/cbmb_allocator.h"

namespace espreso {
namespace gpu {
namespace mgm {

    struct _device {};

    struct _queue {};

    device get_device_by_mpi(int mpi_rank, int mpi_size) { return device{}; }

    void init_gpu(device & d) {}

    void set_device(device & d) {}

    void queue_create(queue & q) {}

    void queue_destroy(queue & q) {}

    void queue_async_barrier(const std::vector<queue> & waitfor, const std::vector<queue> & waitin) {}

    void queue_wait(queue & q) {}

    void device_wait() {}

    size_t get_device_memory_capacity() { return 0; }

    void * memalloc_device(size_t num_bytes) { return nullptr; }

    void memfree_device(void * ptr) {}

    void memalloc_device_max(void * & memory, size_t & memory_size_B, size_t max_needed) {}

    void submit_host_function(queue & q, const std::function<void(void)> & f) {}

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



    #define INSTANTIATE_T_I_AHOST_ADEVICE(T,I,Ahost,Adevice) \
    template void copy_submit_h2d<T,I,Adevice,Ahost  >(queue & q, Vector_Dense<T,I,Adevice> & output, const Vector_Dense<T,I,Ahost>   & input); \
    template void copy_submit_d2h<T,I,Ahost,  Adevice>(queue & q, Vector_Dense<T,I,Ahost>   & output, const Vector_Dense<T,I,Adevice> & input); \
    template void copy_submit_h2d<T,I,Adevice,Ahost  >(queue & q, Matrix_Dense<T,I,Adevice> & output, const Matrix_Dense<T,I,Ahost>   & input); \
    template void copy_submit_d2h<T,I,Ahost,  Adevice>(queue & q, Matrix_Dense<T,I,Ahost>   & output, const Matrix_Dense<T,I,Adevice> & input); \
    template void copy_submit_h2d<T,I,Adevice,Ahost  >(queue & q, Matrix_CSR<T,I,Adevice> & output, const Matrix_CSR<T,I,Ahost>   & input, bool copy_pattern = true, bool copy_vals = true); \
    template void copy_submit_d2h<T,I,Ahost,  Adevice>(queue & q, Matrix_CSR<T,I,Ahost>   & output, const Matrix_CSR<T,I,Adevice> & input, bool copy_pattern = true, bool copy_vals = true);

    #define INSTANTIATE_T_I(T,I) \
    template void copy_submit_h2d<T,I>(queue & q, T * dst, const T * src, I num_elements); \
    template void copy_submit_d2h<T,I>(queue & q, T * dst, const T * src, I num_elements); \
    INSTANTIATE_T_I_AHOST_ADEVICE(T, I, mgm::Ah,       mgm::Ad) \
    INSTANTIATE_T_I_AHOST_ADEVICE(T, I, mgm::Ah,       cbmba_d) \
    INSTANTIATE_T_I_AHOST_ADEVICE(T, I, cpu_allocator, mgm::Ad) \
    INSTANTIATE_T_I_AHOST_ADEVICE(T, I, cpu_allocator, cbmba_d)

    #define INSTANTIATE_T(T) \
    INSTANTIATE_T_I(T, int32_t) \
    /* INSTANTIATE_T_I(T, int64_t) */

    // INSTANTIATE_T(float)
    INSTANTIATE_T(double)
    // INSTANTIATE_T(std::complex<float>)
    // INSTANTIATE_T(std::complex<double>)
    // INSTANTIATE_T(float*)
    INSTANTIATE_T(double*)
    // INSTANTIATE_T(std::complex<float>*)
    // INSTANTIATE_T(std::complex<double>*)
    INSTANTIATE_T(int32_t)
    INSTANTIATE_T(int32_t*)
    // INSTANTIATE_T(int64_t)
    // INSTANTIATE_T(int64_t*)

    #undef INSTANTIATE_T
    #undef INSTANTIATE_T_I
    #undef INSTANTIATE_T_I_AHOST_ADEVICE

}
}
}

#endif
