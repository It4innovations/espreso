
#ifdef HAVE_ONEAPI
#ifdef ESPRESO_USE_WRAPPER_GPU_ONEAPI

#include "gpu/gpu_management.h"
#include "common_oneapi_mgm.h"
#include "basis/utilities/cbmb_allocator.h"
#pragma clang diagnostic push
#pragma clang diagnostic ignored "-Wall"
#pragma clang diagnostic ignored "-Wdeprecated-declarations"
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

    bool is_available()
    {
        std::vector<sycl::device> all_gpus = sycl::device::get_devices(sycl::info::device_type::gpu);
        std::vector<sycl::device> gpus_levelzero;
        for(sycl::device & gpu : all_gpus) {
            if(gpu.get_backend() == sycl::backend::ext_oneapi_level_zero) {
                gpus_levelzero.push_back(gpu);
            }
        }
        return (gpus_levelzero.size() > 0);
    }

    device get_device_by_mpi(int mpi_rank, int mpi_size)
    {
        std::vector<sycl::device> all_gpus = sycl::device::get_devices(sycl::info::device_type::gpu);
        std::vector<sycl::device> gpus_levelzero;
        for(sycl::device & gpu : all_gpus) {
            if(gpu.get_backend() == sycl::backend::ext_oneapi_level_zero) {
                gpus_levelzero.push_back(gpu);
            }
        }

        int my_gpu_idx = gpu::mgm::_internal_get_local_gpu_idx(device_count);

        sycl::device & selected_dev = gpus_levelzero[my_gpu_idx];
        sycl::device d = sycl::device(selected_dev);
        sycl::context c = sycl::context(d);
        device dev = std::make_shared<_device>(d, c);
        return dev;
    }

    void init_gpu(device & d)
    {
        d->mem_allocated = 0;
    }

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

    size_t get_device_memory_free()
    {
        // this works wrong on 2-stack GPU, dont use
        // return default_device->d.get_info<sycl::ext::intel::info::device::free_memory>();

        size_t capacity = get_device_memory_capacity();
        size_t capacity_with_margin = (capacity * 95) / 100;
        size_t allocated = default_device->mem_allocated;
        size_t allocated_with_margin = (allocated * 110) / 100;
        if(capacity_with_margin > allocated_with_margin) {
            return capacity_with_margin - allocated_with_margin;
        }
        else {
            return 0;
        }
    }

    size_t get_natural_pitch_align()
    {
        return 512;
    }

    void * memalloc_device(size_t num_bytes)
    {
        if(num_bytes == 0) return nullptr;

        std::lock_guard<std::mutex> lock(default_device->mtx_alloc);

        void * ptr = sycl::malloc_device(num_bytes, default_device->d, default_device->c);
        if(ptr == nullptr) {
            eslog::error("memalloc_device: failed to allocate memory\n");
        }
        default_device->mem_allocated += num_bytes;
        default_device->alloc_sizes.insert({ptr,num_bytes});
        return ptr;
    }

    void * memalloc_device_2d(size_t num_chunks, size_t bytes_per_chunk, size_t & pitch)
    {
        constexpr size_t align = 512;
        pitch = ((bytes_per_chunk - 1) / align + 1) * align;
        size_t total_size = num_chunks * pitch;
        return memalloc_device(total_size);
    }

    void memfree_device(void * ptr)
    {
        if(ptr == nullptr) return;

        std::lock_guard<std::mutex> lock(default_device->mtx_alloc);

        auto it = default_device->alloc_sizes.find(ptr);
        size_t num_bytes = it->second;
        default_device->alloc_sizes.erase(it);
        default_device->mem_allocated -= num_bytes;
        sycl::free(ptr, default_device->c);
    }

    void memalloc_device_max(void * & memory, size_t & memory_size_B, size_t max_needed)
    {
        size_t keep_free_percent = 5;
        size_t size_free = get_device_memory_free();
        size_t can_allocate_max = ((100 - keep_free_percent) * size_free) / 100;
        memory_size_B = std::min(can_allocate_max, max_needed);
        memory = memalloc_device(memory_size_B);
    }

    void * memalloc_hostpinned(size_t num_bytes)
    {
        if(num_bytes == 0) return nullptr;

        // void * ptr = sycl::malloc_host(num_bytes, default_device->c);
        // if(ptr == nullptr) {
        //     eslog::error("memalloc_hostpinned: failed to allocate memory\n");
        // }

        void * ptr = malloc(num_bytes);
        sycl::ext::oneapi::experimental::prepare_for_device_copy(ptr, num_bytes, default_device->c);

        return ptr;
    }

    void memfree_hostpinned(void * ptr)
    {
        // sycl::free(ptr, default_device->c);

        sycl::ext::oneapi::experimental::release_from_device_copy(ptr, default_device->c);
        free(ptr);
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

    template<typename T>
    void copy_submit(queue & q, PermutationView_new<T> & src, PermutationView_new<T> & dst)
    {
        if(src.size != dst.size) eslog::error("copy submit: output permutation has wrong dimensions\n");
        copy_submit(q, dst.src_to_dst, src.src_to_dst, src.size);
        copy_submit(q, dst.dst_to_src, src.dst_to_src, src.size);
    }

    template<typename T>
    void copy_submit(queue & q, VectorDenseView_new<T> & src, VectorDenseView_new<T> & dst)
    {
        if(src.size != dst.size) eslog::error("copy submit: incompatible sizes\n");

        q->q.template copy<T>(src.vals, dst.vals, dst.size);
    }

    template<typename T, typename I>
    void copy_submit(queue & q, MultiVectorDenseView_new<T,I> & src, MultiVectorDenseView_new<T,I> & dst, bool copy_pattern, bool copy_values)
    {
        if(src.num_vectors != dst.num_vectors || src.size != dst.size) eslog::error("copy submit: incompatible dimensions\n");

        if(copy_pattern) q->q.template copy<I>(src.offsets, dst.offsets, src.num_vectors + 1);
        if(copy_values) q->q.template copy<T>(src.vals, dst.vals, src.size);
    }

    template<typename T>
    void copy_submit(queue & q, MatrixDenseView_new<T> & src, MatrixDenseView_new<T> & dst)
    {
        if(src.order != dst.order) eslog::error("copy submit: orders must match\n");
        if(src.nrows != dst.nrows || src.ncols != dst.ncols) eslog::error("copy submit: incompatible matrix dimensions\n");

        dpct::async_dpct_memcpy(dst.vals, dst.ld * sizeof(T), src.vals, src.ld * sizeof(T), dst.get_size_secdary() * sizeof(T), dst.get_size_primary(), dpct::memcpy_direction::automatic, q->q);
    }

    template<typename T, typename I>
    void copy_submit(queue & q, MatrixCsxView_new<T,I> & src, MatrixCsxView_new<T,I> & dst, bool copy_pattern, bool copy_vals)
    {
        if(src.order != dst.order) eslog::error("copy submit: orders must match\n");
        if(src.nrows != dst.nrows || src.ncols != dst.ncols || src.nnz != dst.nnz) eslog::error("copy submit: incompatible matrix dimensions\n");

        if(copy_pattern) q->q.template copy<I>(src.ptrs, dst.ptrs, dst.get_size_primary() + 1);
        if(copy_pattern) q->q.template copy<I>(src.idxs, dst.idxs, dst.nnz);
        if(copy_vals) q->q.template copy<T>(src.vals, dst.vals, dst.nnz);
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
#endif
