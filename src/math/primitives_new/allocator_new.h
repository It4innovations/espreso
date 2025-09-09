
#ifndef SRC_MATH_PRIMITIVES_NEW_ALLOCATOR_NEW_H_
#define SRC_MATH_PRIMITIVES_NEW_ALLOCATOR_NEW_H_

#include "gpu/gpu_management.h"
#include "basis/utilities/cbmb_allocator.h"



namespace espreso {



class AllocatorCPU_new : public Allocator_new
{
private:
    static constexpr size_t align = 64; // cache line alignment
    static AllocatorCPU_new singleton;
public:
    AllocatorCPU_new() {}
    virtual ~AllocatorCPU_new() {}
    virtual void * alloc(size_t num_bytes) override
    {
        return aligned_alloc(align, num_bytes);
    }
    virtual void free(void * & ptr) override
    {
        ::free(ptr);
        ptr = nullptr;
    }
    virtual bool is_data_accessible_cpu() override
    {
        return true;
    }
    virtual bool is_data_accessible_gpu() override
    {
        return false;
    }
    virtual size_t get_align() override
    {
        return align;
    }
public:
    static AllocatorCPU_new * get_singleton()
    {
        return &AllocatorCPU_new::singleton;
    }
};



class AllocatorGPU_new : public Allocator_new
{
private:
    static AllocatorGPU_new singleton;
public:
    AllocatorGPU_new() {}
    virtual ~AllocatorGPU_new() {}
    virtual void * alloc(size_t num_bytes) override
    {
        return gpu::mgm::memalloc_device(num_bytes);
    }
    virtual void free(void * & ptr) override
    {
        gpu::mgm::memfree_device(ptr);
        ptr = nullptr;
    }
    virtual bool is_data_accessible_cpu() override
    {
        return false;
    }
    virtual bool is_data_accessible_gpu() override
    {
        return true;
    }
    virtual size_t get_align() override
    {
        return gpu::mgm::get_natural_pitch_align();
    }
public:
    static AllocatorGPU_new * get_singleton()
    {
        return &AllocatorGPU_new::singleton;
    }
};



class AllocatorHostPinned_new : public Allocator_new
{
private:
    static AllocatorHostPinned_new singleton;
public:
    AllocatorHostPinned_new() {}
    virtual ~AllocatorHostPinned_new() {}
    virtual void * alloc(size_t num_bytes) override
    {
        return gpu::mgm::memalloc_hostpinned(num_bytes);
    }
    virtual void free(void * & ptr) override
    {
        gpu::mgm::memfree_hostpinned(ptr);
        ptr = nullptr;
    }
    virtual bool is_data_accessible_cpu() override
    {
        return true;
    }
    virtual bool is_data_accessible_gpu() override
    {
        return false;
    }
    virtual size_t get_align() override
    {
        return 64;
    }
public:
    static AllocatorHostPinned_new * get_singleton()
    {
        return &AllocatorHostPinned_new::singleton;
    }
};



class AllocatorArena_new : public Allocator_new
{
private:
    Allocator_new * origin_ator;
    char * start_ptr = nullptr;
    size_t capacity = 0;
    size_t curr_used = 0;
    size_t align_B;
public:
    AllocatorArena_new(Allocator_new * origin_ator_) : origin_ator(origin_ator_)
    {
        align_B = origin_ator->get_align();
    }
    virtual ~AllocatorArena_new() {}
    void set(void * ptr_, size_t capacity_)
    {
        if(start_ptr != nullptr) eslog::error("arena allocator has already been set\n");
        if((uintptr_t)ptr_ % align_B != 0) eslog::error("arena buffer pointer must be a multiple of align\n");

        start_ptr = (char*)ptr_;
        capacity = capacity_;
        curr_used = 0;
    }
    void clear()
    {
        curr_used = 0;
    }
    void unset()
    {
        start_ptr = nullptr;
    }
    virtual void * alloc(size_t num_bytes) override
    {
        if(num_bytes == 0) return nullptr;
        if(start_ptr == nullptr) eslog::error("arena allocator has not been set yet\n");

        char * ptr = nullptr;
        #pragma omp critical(espreso_AllocatorArena_new_alloc)
        {
            ptr = start_ptr + curr_used;
            curr_used += num_bytes;
            if(curr_used > capacity) {
                eslog::error("arena allocator exceeded its capacity\n");
            }
            curr_used = utils::round_up(curr_used, align_B);
        }
        return ptr;
    }
    virtual void free(void * & ptr) override
    {
        ptr = nullptr;
    }
    virtual bool is_data_accessible_cpu() override
    {
        return origin_ator->is_data_accessible_cpu();
    }
    virtual bool is_data_accessible_gpu() override
    {
        return origin_ator->is_data_accessible_gpu();
    }
    virtual size_t get_align() override
    {
        return align_B;
    }
    size_t get_remaining_capacity()
    {
        return utils::round_down(capacity - curr_used, align_B);
    }
};



class AllocatorSinglePointer_new : public Allocator_new
{
private:
    Allocator_new * origin_ator;
    void * pointer = nullptr;
    size_t capacity = 0;
    size_t align_B;
public:
    AllocatorSinglePointer_new(Allocator_new * origin_ator_) : origin_ator(origin_ator_)
    {
        align_B = origin_ator->get_align();
    }
    virtual ~AllocatorSinglePointer_new() {}
    void set(void * ptr_, size_t capacity_)
    {
        if(pointer != nullptr) eslog::error("singlepointer allocator has already been set\n");
        if((uintptr_t)ptr_ % align_B != 0) eslog::error("singlepointer buffer pointer must be a multiple of align\n");

        pointer = ptr_;
        capacity = capacity_;
    }
    void unset()
    {
        pointer = nullptr;
    }
    virtual void * alloc(size_t num_bytes) override
    {
        if(num_bytes == 0) return nullptr;
        if(pointer == nullptr) eslog::error("singlepointer allocator has not been set yet\n");
        if(num_bytes > capacity) {
            eslog::error("singlepointer allocator exceeded its capacity (%zu > %zu)\n", num_bytes, capacity);
        }
        return pointer;
    }
    virtual void free(void * & ptr) override
    {
        ptr = nullptr;
    }
    virtual bool is_data_accessible_cpu() override
    {
        return origin_ator->is_data_accessible_cpu();
    }
    virtual bool is_data_accessible_gpu() override
    {
        return origin_ator->is_data_accessible_gpu();
    }
    virtual size_t get_align() override
    {
        return align_B;
    }
};



class AllocatorCBMB_new : public Allocator_new
{
private:
    Allocator_new * origin_ator;
    size_t align_B;
public:
    cbmba_resource resource;
public:
    AllocatorCBMB_new(Allocator_new * origin_ator_, void * memory, size_t capacity) :  origin_ator(origin_ator_), resource(memory, capacity)
    {
        align_B = origin_ator->get_align();
    }
    AllocatorCBMB_new(const AllocatorCBMB_new & other) = delete;
    AllocatorCBMB_new(AllocatorCBMB_new && other) = delete;
    AllocatorCBMB_new & operator=(const AllocatorCBMB_new & other) = delete;
    AllocatorCBMB_new & operator=(AllocatorCBMB_new && other) = delete;
    virtual ~AllocatorCBMB_new() {}
    virtual void * alloc(size_t num_bytes) override
    {
        return resource.allocate(num_bytes, align_B);
    }
    virtual void free(void * & ptr) override
    {
        resource.deallocate(ptr);
        ptr = nullptr;
    }
    virtual bool is_data_accessible_cpu() override
    {
        return origin_ator->is_data_accessible_cpu();
    }
    virtual bool is_data_accessible_gpu() override
    {
        return origin_ator->is_data_accessible_gpu();
    }
    virtual size_t get_align() override
    {
        return align_B;
    }
};



}



#endif /* SRC_MATH_PRIMITIVES_NEW_ALLOCATOR_NEW_H_ */
