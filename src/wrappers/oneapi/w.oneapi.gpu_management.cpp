
#ifdef HAVE_ONEAPI

#include "gpu/gpu_management.h"
#include "w.oneapi.gpu_management.h"
#include "basis/utilities/cbmb_allocator.h"
#pragma clang diagnostic push
#pragma clang diagnostic ignored "-Wunused-function"
#pragma clang diagnostic ignored "-Wunneeded-internal-declaration"
#pragma clang diagnostic ignored "-Wunused-parameter"
#pragma clang diagnostic ignored "-Wsign-compare"
#include <dpct/dpct.hpp>
#pragma clang diagnostic pop

namespace espreso {
namespace gpu {
namespace mgm {

    device default_device = device{};

    gpu_wrapper_impl get_implementation()
    {
        return gpu_wrapper_impl::ONEAPI;
    }

    device get_device_by_mpi(int /*mpi_rank*/, int /*mpi_size*/)
    {
        // TODO: select with knowledge of mpi
        sycl::device d = sycl::device(sycl::gpu_selector_v);
        sycl::context c = sycl::context(d);
        device dev = std::make_shared<_device>(d, c);
        return dev;
    }

    void init_gpu(device & /*d*/) {}

    void set_device(device & d)
    {
        default_device = d;
    }

    void queue_create(queue & q)
    {
        q = std::make_shared<_queue>(sycl::queue(default_device->c, default_device->d, sycl::property::queue::in_order()));
        q->d = default_device;
        q->d->qs.push_back(q);
    }

    void queue_destroy(queue & q)
    {
        std::vector<queue> & dqs = q->d->qs;
        dqs.erase(std::find(dqs.begin(), dqs.end(), q));
        q.reset();
    }

    void queue_async_barrier(const std::vector<queue> & waitfor, const std::vector<queue> & waitin)
    {
        std::vector<sycl::event> events;
        events.reserve(waitfor.size());
        for(const queue & q : waitfor) {
            events.push_back(q->q.single_task([](){}));
        }
        for(const queue & q : waitin) {
            q->q.single_task(events, [](){});
        }
    }

    void queue_wait(queue & q)
    {
        q->q.wait();
    }

    void device_wait()
    {
        for(queue & q : default_device->qs) {
            queue_wait(q);
        }
    }

    size_t get_device_memory_capacity()
    {
        size_t capacity = default_device->d.get_info<sycl::info::device::global_mem_size>();
        return capacity;
    }

    void * memalloc_device(size_t num_bytes)
    {
        void * ptr = sycl::malloc_device(num_bytes, default_device->d, default_device->c);
        return ptr;
    }

    void memfree_device(void * ptr)
    {
        sycl::free(ptr, default_device->c);
    }

    void memalloc_device_max(void * & memory, size_t & memory_size_B, size_t max_needed)
    {
        size_t keep_free_percent = 5;
        size_t size_free = default_device->d.get_info<sycl::ext::intel::info::device::free_memory>();
        size_t can_allocate_max = ((100 - keep_free_percent) * size_free) / 100;
        memory_size_B = std::min(can_allocate_max, max_needed);
        memory = memalloc_device(memory_size_B);
    }

    void * memalloc_hostpinned(size_t num_bytes)
    {
        void * ptr = sycl::malloc_host(num_bytes, default_device->c);
        return ptr;
    }

    void memfree_hostpinned(void * ptr)
    {
        sycl::free(ptr, default_device->c);
    }

    void submit_host_function(queue & q, const std::function<void(void)> & f)
    {
        q->q.submit([&](sycl::handler & cgh){
            cgh.host_task([=](){
                f();
            });
        });
    }

    template<typename T>
    void copy_submit(queue & q, T * dst, T const * src, size_t num_elements)
    {
        q->q.copy<T>(src, dst, num_elements);
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
        if(output.get_ld() == input.get_ld()) {
            copy_submit(q, output.vals, input.vals, input.nrows * input.get_ld());
        }
        else {
            dpct::async_dpct_memcpy(output.vals, output.get_ld() * sizeof(T), input.vals, input.get_ld() * sizeof(T), input.ncols * sizeof(T), input.nrows, dpct::memcpy_direction::automatic, q->q);
        }
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
        q->q.fill(ptr, val, num_bytes);
    }

}
}
}

#include "gpu/gpu_management.inst.hpp"

#endif
