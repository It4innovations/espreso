
#ifndef SRC_MATH_PRIMITIVES_NEW_ALLOCATOR_NEW_BASE_H_
#define SRC_MATH_PRIMITIVES_NEW_ALLOCATOR_NEW_BASE_H_

#include <cstddef>
#include "esinfo/eslog.h"



namespace espreso {



// TODO idea about cpu-gpu location
// have separate "ability to access data" and "efficient access to data"
// e.g. GPU can access cudaMallocHost-allocated data, but it is not efficient
// also I should solve multi-gpu environments

class Allocator_new
{
public:
    Allocator_new() {}
    virtual ~Allocator_new() {}
    virtual void * alloc(size_t num_bytes) = 0;
    virtual void free(void * & ptr) = 0; // must set ptr to nullptr
    virtual bool is_data_accessible_cpu() = 0;
    virtual bool is_data_accessible_gpu() = 0;
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
    virtual bool is_data_accessible_cpu() override
    {
        return on_cpu;
    }
    virtual bool is_data_accessible_gpu() override
    {
        return on_gpu;
    }
public:
    static AllocatorDummy_new * get_singleton(bool cpu, bool gpu)
    {
        if(!cpu && !gpu) return &singleton_ff;
        if(!cpu &&  gpu) return &singleton_ft;
        if( cpu && !gpu) return &singleton_tf;
        if( cpu &&  gpu) return &singleton_tt;
        eslog::error("unreachable code\n");
    }
};



}



#endif /* SRC_MATH_PRIMITIVES_NEW_ALLOCATOR_NEW_BASE_H_ */
