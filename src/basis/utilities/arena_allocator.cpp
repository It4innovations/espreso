
#include "arena_allocator.h"

#include "esinfo/eslog.hpp"



namespace espreso {

    arena_allocator_resource::~arena_allocator_resource() = default;

    arena_allocator_resource::arena_allocator_resource(void * mem_, size_t capacity_, size_t align_B_) : mem((char*)mem_), capacity(capacity_), align_B(align_B_)
    {
        if(reinterpret_cast<std::uintptr_t>(mem) % align_B != 0) {
            eslog::error("arena_allocator: wrong pointer alignment\n");
        }
    }

    void * arena_allocator_resource::allocate(size_t num_bytes)
    {
        if(num_bytes == 0) return nullptr;

        std::lock_guard lock(mtx);

        void * ptr = mem + curr_offset;
        curr_offset += num_bytes;
        if(curr_offset > capacity) {
            eslog::error("arena_allocator: capacity is full\n");
        }
        curr_offset = ((curr_offset + 1) / align_B + 1) * align_B;

        num_allocations++;
        return ptr;
    }

    void arena_allocator_resource::deallocate(void * & ptr)
    {
        if(ptr == nullptr) return;

        std::lock_guard lock(mtx);

        num_allocations--;
        ptr = nullptr;
    }

    void arena_allocator_resource::reset_arena()
    {
        std::lock_guard lock(mtx);

        if(num_allocations != 0) {
            eslog::error("arena_allocator: mismatch between allocations and deallocations, %zu allocations not deallocated\n", num_allocations);
        }

        num_allocations = 0;
        curr_offset = 0;
    }

    size_t arena_allocator_resource::get_free_size()
    {
        return capacity - curr_offset;
    }

    size_t arena_allocator_resource::get_capacity()
    {
        return capacity;
    }

}
