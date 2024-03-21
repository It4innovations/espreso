
#ifndef SRC_BASIS_UTILITIES_CBMB_ALLOCATOR_H_
#define SRC_BASIS_UTILITIES_CBMB_ALLOCATOR_H_

#include <algorithm>
#include <list>
#include <mutex>
#include <condition_variable>
#include <omp.h>

#include "esinfo/eslog.hpp"



// Cyclic Buffer Mulithreaded Blocking Allocator (and deallocator)
// not a real C++-conforming allocator, but just my little implementation
// kind of assumes that the memory will be freed approximately in the same order as it was allocated
// the memory pool is owned by the user, but used by the cbmba and should not be touched during the lifetime of the cbmba_resource object

// for a given allocate-deallocate region, all allocations should be in a critical region to avoid deadlocks
//     allocate the memory inside cbmba_resource::do_transaction([](){ /* lambda */ }), all will be managed for you
// beware of fragmentation - there can be an allocated block of memory in the middle the pool, which leaves not enough space on the sides for a new larger allocation
//     to guarantee successfull allocations (in this problem scenario), the total memory that is required for a single loop iteration should be at most half of the pool size. I think. Or the largest chunk of memory must be smaller than half of remaining memory
// the allocated memory should be freed asap to allow it to be used in other memory requests. Allocating a chunk of memory for the lifetime of the object will deadlock, since I use the buffer in a circular way -- unless the memory block is freed, no allocations located after it can occur



namespace espreso {

    class cbmba_resource
    {
    private:
        struct mem_block
        {
            size_t start;
            size_t start_aligned;
            mem_block(size_t s, size_t sa) : start(s), start_aligned(sa) {}
        };
        std::mutex mtx_transaction;
        std::mutex mtx_blocks;
        std::condition_variable cv;
        char * memory_pool;
        size_t memory_pool_size;
        size_t volatile start_empty;
        size_t volatile start_full;
        std::list<mem_block> full_blocks;
        size_t transaction_end_fail;
    public:
        cbmba_resource(void * memory_pool_, size_t memory_pool_size_)
            : memory_pool(reinterpret_cast<char*>(memory_pool_))
            , memory_pool_size(memory_pool_size_)
            , start_empty(0)
            , start_full(0)
            , transaction_end_fail(~(size_t)0)
        {
        }
        ~cbmba_resource()
        {
        }
        cbmba_resource() = delete;
        cbmba_resource(cbmba_resource const &) = delete;
        cbmba_resource(cbmba_resource &&) = delete;
        cbmba_resource & operator=(cbmba_resource const &) = delete;
        cbmba_resource & operator=(cbmba_resource &&) = delete;
        size_t get_capacity() { return memory_pool_size; }
        void * allocate(size_t num_bytes, size_t align = 1)
        {
            if(num_bytes == 0) num_bytes = 1;
            if(num_bytes > memory_pool_size) eslog::error("Not enough memory in the pool, capacity is %zu B = %zu MiB\n", memory_pool_size, memory_pool_size >> 20);
            if(memory_pool_size % align != 0) eslog::error("Align has to divide memory_pool_size\n");
            if(reinterpret_cast<std::uintptr_t>(memory_pool) % align != 0) eslog::error("Memory pool has to be aligned the same way as is the requested align\n");

            size_t new_block_start;
            size_t new_block_start_aligned;
            size_t new_block_end;
            {
                std::unique_lock<std::mutex> lk(mtx_blocks);

                new_block_start = start_empty;
                new_block_start_aligned = ((new_block_start - 1) / align + 1) * align;
                new_block_end = new_block_start_aligned + num_bytes;
                if((new_block_end - 1) / memory_pool_size > new_block_start / memory_pool_size) // if the block would cross the pool size boundary
                {
                    new_block_start = ((start_empty - 1) / memory_pool_size + 1) * memory_pool_size;
                    new_block_start_aligned = ((new_block_start - 1) / align + 1) * align;
                    new_block_end = new_block_start_aligned + num_bytes;
                }

                full_blocks.emplace_back(new_block_start, new_block_start_aligned);
                start_empty = new_block_end;

                if(start_empty > transaction_end_fail) eslog::error("More memory requested in a single transaction then what is available in the pool. Capacity is %zu B = %zu MiB\n", memory_pool_size, memory_pool_size >> 20);

                // printf("Waiting for memory, %zu B, this=%p, new_block_end=%zu, sf+mps=%zu\n", num_bytes, this, new_block_end, start_full + memory_pool_size);
                cv.wait(lk, [&]{ return new_block_end <= start_full + memory_pool_size; });
                // printf("Memory granted,     %zu B, this=%p, new_block_end=%zu, sf+mps=%zu\n", num_bytes, this, new_block_end, start_full + memory_pool_size);
            }

            return memory_pool + new_block_start_aligned % memory_pool_size;
        }
        void deallocate(void * ptr)
        {
            if(ptr == nullptr) return;
            size_t mem_block_start_aligned = reinterpret_cast<char*>(ptr) - memory_pool;

            {
                std::lock_guard<std::mutex> lk(mtx_blocks);

                // printf("Deallocating, this=%p, ptr=%p, pos=%zu\n", this, ptr, mem_block_start_aligned);
                // remove_from_full_blocks_vector(mem_block_start_aligned);
                auto it = std::find_if(full_blocks.begin(), full_blocks.end(), [=](mem_block const & mb){return mb.start_aligned % memory_pool_size == mem_block_start_aligned;});
                if(it == full_blocks.end()) eslog::error("Bad remove from full blocks vector\n");
                full_blocks.erase(it);

                if(full_blocks.size() == 0) start_full = start_empty = 0;
                else start_full = full_blocks.front().start;
            }
            cv.notify_all();
        }
        void print_full_blocks()
        {
            printf("cbmba_resource full blocks:");
            for(mem_block & mb : full_blocks) printf(" %zu/%zu", mb.start, mb.start % memory_pool_size);
            printf("\n");
        }
        void do_transaction(const std::function<void(void)> & f)
        {
            std::lock_guard<std::mutex> lk(mtx_transaction);
            // c must have an operator() with no parameters and must only do allocations using this resource
            transaction_end_fail = start_empty + memory_pool_size;
            f();
            transaction_end_fail = ~(size_t)0;
        }
    };



    template<bool ad, bool ah>
    class cbmba
    {
    public:
        static constexpr bool is_data_device_accessible = ad;
        static constexpr bool is_data_host_accessible = ah;
        static constexpr bool always_equal = true; // for my usage of this class it holds, but not in general
    private:
        cbmba_resource & resource;
        size_t align_B;
    public:
        cbmba() = delete;
        cbmba(cbmba<ad,ah> const &) = default;
        cbmba(cbmba<ad,ah> &&) = default;
        ~cbmba() = default;
        cbmba<ad,ah> & operator=(cbmba<ad,ah> const &) = default;
        cbmba<ad,ah> & operator=(cbmba<ad,ah> &&) = default;
    public:
        cbmba(cbmba_resource & res, size_t align_B_) : resource(res), align_B(align_B_) {}
        void * allocate(size_t num_bytes)
        {
            return resource.allocate(num_bytes, align_B);
        }
        template<typename T>
        T * allocate(size_t count)
        {
            return reinterpret_cast<T*>(allocate(count * sizeof(T)));
        }
        template<typename T>
        void deallocate(T * ptr)
        {
            resource.deallocate(ptr);
        }
    };

    using cbmba_d = cbmba<true,false>;
    using cbmba_h = cbmba<false,true>;

}

#endif /* SRC_BASIS_UTILITIES_CBMB_ALLOCATOR_H_ */
