
#ifdef HAVE_ONEAPI

#include "gpu/gpu_management.h"
#include "w.oneapi.gpu_management.h"
#include "basis/utilities/cbmb_allocator.h"

namespace espreso {
namespace gpu {
namespace mgm {

    device default_device = device{};

    gpu_wrapper_impl get_implementation()
    {
        return gpu_wrapper_impl::ONEAPI;
    }

    device get_device_by_mpi(int mpi_rank, int mpi_size)
    {
        // TODO: select with knowledge of mpi
        device d = std::make_shared<_device>(sycl::device(sycl::gpu_selector_v));
        return d;
    }

    void init_gpu(device & d) {}

    void set_device(device & d)
    {
        default_device = d;
    }

    void queue_create(queue & q)
    {
        q = std::make_shared<_queue>(sycl::queue(default_device->d, sycl::property::queue::in_order()));
        default_device->qs.push_back(q);
    }

    void queue_destroy(queue & q)
    {
        std::vector<queue> & dqs = q->d->qs;
        dqs.erase(std::find(dqs.begin(), dqs.end(), q));
        q.reset();
    }

    void queue_async_barrier(const std::vector<queue> & waitfor, const std::vector<queue> & waitin)
    {
        throw std::runtime_error("unimplemented");
    }

    void queue_wait(queue & q)
    {
        q->q.wait();
    }

    void device_wait()
    {
        for(queue & q : default_device->qs)
        {
            queue_wait(q);
        }
    }

    void event_create(event & e) {}

    void event_destroy(event & e) {}

    void event_record(event & e, queue & q) {}

    void event_wait(event & e) {}

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

#include "gpu/gpu_management.inst.hpp"

#endif
