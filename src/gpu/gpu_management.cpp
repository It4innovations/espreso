
#ifndef ESPRESO_USE_WRAPPER_GPU_CUDA
#ifndef ESPRESO_USE_WRAPPER_GPU_ROCM
#ifndef ESPRESO_USE_WRAPPER_GPU_ONEAPI

#include "gpu_management.h"
#include "basis/utilities/cbmb_allocator.h"

namespace espreso {
namespace gpu {
namespace mgm {

    gpu_wrapper_impl get_implementation()
    {
        return gpu_wrapper_impl::NONE;
    }

    bool is_available() { return false; }

    struct _device {};

    struct _queue {};

    device get_device_by_mpi(int /*mpi_rank*/, int /*mpi_size*/) { eslog::error("calling empty GPU wrapper.\n"); }

    void init_gpu(device & /*d*/) { eslog::error("calling empty GPU wrapper.\n"); }

    void set_device(device & /*d*/) { eslog::error("calling empty GPU wrapper.\n"); }

    void queue_create(queue & /*q*/) { eslog::error("calling empty GPU wrapper.\n"); }

    void queue_destroy(queue & /*q*/) {}

    void queue_async_barrier(const std::vector<queue> & /*waitfor*/, const std::vector<queue> & /*waitin*/) {}

    void queue_wait(queue & /*q*/) {}

    void device_wait() {}

    size_t get_device_memory_capacity() { return 0; }

    size_t get_device_memory_free() { return 0; }

    size_t get_natural_pitch_align() { return 1; }

    void * memalloc_device(size_t /*num_bytes*/) { return nullptr; }

    void * memalloc_device_2d(size_t /*num_chunks*/, size_t /*bytes_per_chunk*/, size_t & /*pitch*/) { return nullptr; }

    void memfree_device(void * /*ptr*/) {}

    void memalloc_device_max(void * & /*memory*/, size_t & /*memory_size_B*/, size_t /*max_needed*/) {}

    void * memalloc_hostpinned(size_t /*num_bytes*/) { return nullptr; }

    void memfree_hostpinned(void * /*ptr*/) {}

    void submit_host_function(queue & /*q*/, const std::function<void(void)> & /*f*/) {}

    template<typename T>
    void copy_submit(queue & /*q*/, T * /*dst*/, T const * /*src*/, size_t /*num_elements*/) {}

    template<typename T, typename I, typename Ao, typename Ai>
    void copy_submit(queue & /*q*/, Vector_Dense<T,I,Ao> & /*output*/, const Vector_Dense<T,I,Ai> & /*input*/) {}

    template<typename T, typename I, typename Ao, typename Ai>
    void copy_submit(queue & /*q*/, Matrix_Dense<T,I,Ao> & /*output*/, const Matrix_Dense<T,I,Ai> & /*input*/) {}
    
    template<typename T, typename I, typename Ao, typename Ai>
    void copy_submit(queue & /*q*/, Matrix_CSR<T,I,Ao> & /*output*/, const Matrix_CSR<T,I,Ai> & /*input*/, bool /*copy_pattern*/, bool /*copy_vals*/) {}

    template<typename T>
    void copy_submit(queue & /*q*/, PermutationView_new<T> & /*src*/, PermutationView_new<T> & /*dst*/) {}

    template<typename T>
    void copy_submit(queue & /*q*/, VectorDenseView_new<T> & /*src*/, VectorDenseView_new<T> & /*dst*/) {}

    template<typename T, typename I>
    void copy_submit(queue & /*q*/, MultiVectorDenseView_new<T,I> & /*src*/, MultiVectorDenseView_new<T,I> & /*dst*/, bool /*copy_pattern*/, bool /*copy_vals*/) {}

    template<typename T>
    void copy_submit(queue & /*q*/, MatrixDenseView_new<T> & /*src*/, MatrixDenseView_new<T> & /*dst*/) {}

    template<typename T, typename I>
    void copy_submit(queue & /*q*/, MatrixCsxView_new<T,I> & /*src*/, MatrixCsxView_new<T,I> & /*dst*/, bool /*copy_pattern*/, bool /*copy_vals*/) {}

    void memset_submit(queue & /*q*/, void * /*ptr*/, size_t /*num_bytes*/, char /*val*/) {}

}
}
}

#include "gpu/gpu_management.inst.hpp"

#endif
#endif
#endif
