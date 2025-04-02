
#ifndef SRC_BASIS_UTILITIES_CBMB_ALLOCATOR_H_
#define SRC_BASIS_UTILITIES_CBMB_ALLOCATOR_H_

#include <algorithm>
#include <list>
#include <mutex>
#include <condition_variable>
#include <functional>
#include <omp.h>



// Cyclic Buffer Mulithreaded Blocking Allocator (and deallocator)
// not a real C++-conforming allocator, but just my little implementation
// kind of assumes that the memory will be freed approximately in the same order as it was allocated
// the memory pool is owned by the user, but used by the cbmba and should not be touched during the lifetime of the cbmba_resource object

// for a given allocate-deallocate region, all allocations should be in a critical region to avoid deadlocks
//     allocate the memory inside cbmba_resource::do_transaction([](){ /* lambda */ }), all will be managed for you
// beware of fragmentation - there can be an allocated block of memory in the middle the pool, which leaves not enough space on the sides for a new larger allocation
//     to guarantee successfull allocations (in this problem scenario), the total memory that is required for a single loop iteration (for a single transaction) should be at most half of the pool size. I think. Or the largest chunk of memory must be smaller than half of remaining memory
//     If I detect that more than half of the resource's memory was allocated in a single transaction, the resource transitions into a state where all transactions start at the start of the memory buffer. This makes the caller wait for longer on allocate(), possibly sacrificing performance, but it will at least work
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
        size_t mem_allocated_in_transaction;
        bool do_reset_to_start_after_transaction = false;
        bool in_transaction = false;
    public:
        cbmba_resource(void * memory_pool_, size_t memory_pool_size_);
        ~cbmba_resource();
        cbmba_resource() = delete;
        cbmba_resource(cbmba_resource const &) = delete;
        cbmba_resource(cbmba_resource &&) = delete;
        cbmba_resource & operator=(cbmba_resource const &) = delete;
        cbmba_resource & operator=(cbmba_resource &&) = delete;
        void * allocate(size_t num_bytes, size_t align = 1);
        void deallocate(void * & ptr);
        void print_full_blocks();
        void do_transaction(const std::function<void(void)> & f);
        size_t get_curr_memory_requested_in_transaction() const;
        size_t get_max_capacity() const;
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
        void deallocate(T * & ptr)
        {
            resource.deallocate((void*&)ptr);
        }
    };

    using cbmba_d = cbmba<true,false>;
    using cbmba_h = cbmba<false,true>;

}

#endif /* SRC_BASIS_UTILITIES_CBMB_ALLOCATOR_H_ */
