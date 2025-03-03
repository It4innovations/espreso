
#ifndef SRC_MATH_PRIMITIVES_NEW_ALLOCATOR_NEW_H_
#define SRC_MATH_PRIMITIVES_NEW_ALLOCATOR_NEW_H_

#include "gpu/gpu_management.h"



namespace espreso {



class Allocator_new
{
public:
    Allocator_new() {}
    virtual ~Allocator_new() {}
    virtual void * alloc(size_t num_bytes) = 0;
    virtual void free(void * & ptr) = 0; // must set ptr to nullptr
    virtual bool is_on_cpu() = 0;
    virtual bool is_on_gpu() = 0;
public:
    virtual size_t get_align()
    {
        return 1;
    }
    template<typename T>
    void free(T * & ptr)
    {
        free((void*&)ptr);
    }
    virtual void * alloc_2d(size_t num_chunks, size_t bytes_per_chunk, size_t & pitch)
    {
        pitch = bytes_per_chunk;
        return alloc(num_chunks * pitch);
    }
    template<typename T>
    T * alloc(size_t num_elements)
    {
        return reinterpret_cast<T*>(alloc(num_elements * sizeof(T)));
    }
};



class AllocatorDummy_new : public Allocator_new
{
private:
    bool on_cpu = false;
    bool on_gpu = false;
    static AllocatorDummy_new singleton_ff;
    static AllocatorDummy_new singleton_ft;
    static AllocatorDummy_new singleton_tf;
    static AllocatorDummy_new singleton_tt;
public:
    AllocatorDummy_new(bool cpu, bool gpu) : on_cpu(cpu), on_gpu(gpu) {}
    virtual ~AllocatorDummy_new() {}
    virtual void * alloc(size_t num_bytes) override
    {
        eslog::error("dummy allocator cannot alloc\n");
    }
    virtual void free(void * & ptr) override
    {
        eslog::error("dummy allocator cannot free\n");
    }
    virtual bool is_on_cpu() override
    {
        return on_cpu;
    }
    virtual bool is_on_gpu() override
    {
        return on_gpu;
    }
public:
    AllocatorDummy_new & get_singleton(bool cpu, bool gpu)
    {
        if(!cpu && !gpu) return singleton_ff;
        if(!cpu &&  gpu) return singleton_ft;
        if( cpu && !gpu) return singleton_tf;
        if( cpu &&  gpu) return singleton_tt;
        eslog::error("unreachable code\n");
    }
};



class AllocatorCPU_new : public Allocator_new
{
private:
    static constexpr size_t align = 64;
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
    virtual bool is_on_cpu() override
    {
        return true;
    }
    virtual bool is_on_gpu() override
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
    virtual bool is_on_cpu() override
    {
        return false;
    }
    virtual bool is_on_gpu() override
    {
        return true;
    }
    virtual size_t get_align() override
    {
        eslog::error("todo\n");
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
    virtual bool is_on_cpu() override
    {
        return true;
    }
    virtual bool is_on_gpu() override
    {
        return false;
    }
    virtual size_t get_align() override
    {
        eslog::error("todo\n");
    }
public:
    static AllocatorHostPinned_new * get_singleton()
    {
        return &AllocatorHostPinned_new::singleton;
    }
};



class AllocatorArena_new : public Allocator_new
{

};



class AllocatorCBMB_new : public Allocator_new
{

};



}



#endif /* SRC_MATH_PRIMITIVES_NEW_ALLOCATOR_NEW_H_ */
