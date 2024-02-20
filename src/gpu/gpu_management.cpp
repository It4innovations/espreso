
#if !defined(HAVE_CUDA) && !defined(HAVE_ROCM)

#include "gpu_management.h"
#include "basis/utilities/cbmb_allocator.h"

namespace espreso {
namespace gpu {
namespace mgm {

    gpu_wrapper_impl get_implementation()
    {
        return gpu_wrapper_impl::NONE;
    }

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

    void * memalloc_hostpinned(size_t num_bytes) { return nullptr; }

    void memfree_hostpinned(void * ptr) {}

    void submit_host_function(queue & q, const std::function<void(void)> & f) {}

    template<typename T>
    void copy_submit_h2d(queue & q, T * dst, T const * src, size_t num_elements) {}

    template<typename T>
    void copy_submit_d2h(queue & q, T * dst, T const * src, size_t num_elements) {}

    template<typename T, typename I, typename Ao, typename Ai>
    void copy_submit_d2h(queue & q, Vector_Dense<T,I,Ao> & output, const Vector_Dense<T,I,Ai> & input) {}

    template<typename T, typename I, typename Ao, typename Ai>
    void copy_submit_h2d(queue & q, Vector_Dense<T,I,Ao> & output, const Vector_Dense<T,I,Ai> & input) {}

    template<typename T, typename I, typename Ao, typename Ai>
    void copy_submit_d2h(queue & q, Matrix_Dense<T,I,Ao> & output, const Matrix_Dense<T,I,Ai> & input) {}

    template<typename T, typename I, typename Ao, typename Ai>
    void copy_submit_h2d(queue & q, Matrix_Dense<T,I,Ao> & output, const Matrix_Dense<T,I,Ai> & input) {}
    
    template<typename T, typename I, typename Ao, typename Ai>
    void copy_submit_d2h(queue & q, Matrix_CSR<T,I,Ao> & output, const Matrix_CSR<T,I,Ai> & input, bool copy_pattern, bool copy_vals) {}

    template<typename T, typename I, typename Ao, typename Ai>
    void copy_submit_h2d(queue & q, Matrix_CSR<T,I,Ao> & output, const Matrix_CSR<T,I,Ai> & input, bool copy_pattern, bool copy_vals) {}

    void memset_submit(queue & q, void * ptr, size_t num_bytes, char val) {}

}
}
}

#include "gpu/gpu_management_inst.hpp"

#endif
