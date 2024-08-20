
#ifdef HAVE_ONEAPI

#include "gpu/gpu_management.h"
#include "w.oneapi.gpu_management.h"
#include "basis/utilities/cbmb_allocator.h"
#include <dpct/dpct.hpp>

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
        sycl::device d = sycl::device(sycl::gpu_selector_v);
        sycl::context c = sycl::context(d);
        device dev = std::make_shared<_device>(d, c);
        return dev;
    }

    void init_gpu(device & d) {}

    void set_device(device & d)
    {
        default_device = d;
    }

    void queue_create(queue & q)
    {
        q = std::make_shared<_queue>(sycl::queue(default_device->c, default_device->d, sycl::property::queue::in_order()));
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

    size_t get_device_memory_capacity()
    {
        device.get_info<sycl::info::device::global_mem_size>();
    }

    void * memalloc_device(size_t num_bytes)
    {
        sycl::malloc_device(num_bytes, default_device->d, default_device->c);
    }

    void memfree_device(void * ptr)
    {
        sycl::free(ptr, default_device->c);
    }

    void memalloc_device_max(void * & memory, size_t & memory_size_B, size_t max_needed) {}

    void * memalloc_hostpinned(size_t num_bytes)
    {
        sycl::malloc_host(num_bytes, default_device->c);
    }

    void memfree_hostpinned(void * ptr)
    {
        sycl::free(ptr, default_device->c);
    }

    void submit_host_function(queue & q, const std::function<void(void)> & f)
    {
        q->q.single_task([=](){
            f();
        })
    }

    template<typename T>
    void copy_submit_h2d(queue & q, T * dst, T const * src, size_t num_elements)
    {
        q->q.copy<T>(src, dst, num_elements);
    }

    template<typename T>
    void copy_submit_d2h(queue & q, T * dst, T const * src, size_t num_elements)
    {
        q->q.copy<T>(src, dst, num_elements);
    }

    template<typename T, typename I, typename Ao, typename Ai>
    void copy_submit_d2h(queue & q, Vector_Dense<T,I,Ao> & output, const Vector_Dense<T,I,Ai> & input)
    {
        static_assert(Ao::is_data_host_accessible, "output vector data has to be host accessible");
        static_assert(Ai::is_data_device_accessible, "input vector data has to be device accessible");
        copy_submit_d2h(q, output.vals, input.vals, input.size);
    }

    template<typename T, typename I, typename Ao, typename Ai>
    void copy_submit_h2d(queue & q, Vector_Dense<T,I,Ao> & output, const Vector_Dense<T,I,Ai> & input)
    {
        static_assert(Ao::is_data_device_accessible, "output vector data has to be device accessible");
        static_assert(Ai::is_data_host_accessible, "input vector data has to be host accessible");
        copy_submit_h2d(q, output.vals, input.vals, input.size);
    }

    template<typename T, typename I, typename Ao, typename Ai>
    void copy_submit_d2h(queue & q, Matrix_Dense<T,I,Ao> & output, const Matrix_Dense<T,I,Ai> & input)
    {
        static_assert(Ao::is_data_host_accessible, "output matrix data has to be host accessible");
        static_assert(Ai::is_data_device_accessible, "input matrix data has to be device accessible");
        if(output.get_ld() == input.get_ld()) {
            copy_submit_d2h(q, output.vals, input.vals, input.nrows * input.get_ld());
        }
        else {
            dpct_memcpy(q->q, output.vals, input.vals, output.get_ld() * sizeof(T), input.get_ld() * sizeof(T), input.ncols, input.nrows);
        }
    }

    template<typename T, typename I, typename Ao, typename Ai>
    void copy_submit_h2d(queue & q, Matrix_Dense<T,I,Ao> & output, const Matrix_Dense<T,I,Ai> & input)
    {
        static_assert(Ao::is_data_device_accessible, "output matrix data has to be device accessible");
        static_assert(Ai::is_data_host_accessible, "input matrix data has to be host accessible");
        if(output.get_ld() == input.get_ld()) {
            copy_submit_h2d(q, output.vals, input.vals, input.nrows * input.get_ld());
        }
        else {
            dpct_memcpy(q->q, output.vals, input.vals, output.get_ld() * sizeof(T), input.get_ld() * sizeof(T), input.ncols, input.nrows);
        }
    }
    
    template<typename T, typename I, typename Ao, typename Ai>
    void copy_submit_d2h(queue & q, Matrix_CSR<T,I,Ao> & output, const Matrix_CSR<T,I,Ai> & input, bool copy_pattern, bool copy_vals)
    {
        static_assert(Ao::is_data_host_accessible, "output matrix data has to be host accessible");
        static_assert(Ai::is_data_device_accessible, "input matrix data has to be device accessible");
        copy_submit_d2h(q, output.rows, input.rows, input.nrows+1);
        copy_submit_d2h(q, output.cols, input.cols, input.nnz);
        copy_submit_d2h(q, output.vals, input.vals, input.nnz);
    }

    template<typename T, typename I, typename Ao, typename Ai>
    void copy_submit_h2d(queue & q, Matrix_CSR<T,I,Ao> & output, const Matrix_CSR<T,I,Ai> & input, bool copy_pattern, bool copy_vals)
    {
        static_assert(Ao::is_data_device_accessible, "output matrix data has to be device accessible");
        static_assert(Ai::is_data_host_accessible, "input matrix data has to be host accessible");
        copy_submit_h2d(q, output.rows, input.rows, input.nrows+1);
        copy_submit_h2d(q, output.cols, input.cols, input.nnz);
        copy_submit_h2d(q, output.vals, input.vals, input.nnz);
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
