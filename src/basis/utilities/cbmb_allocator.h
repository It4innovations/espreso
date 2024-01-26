
#ifndef SRC_BASIS_UTILITIES_CBMB_ALLOCATOR_H_
#define SRC_BASIS_UTILITIES_CBMB_ALLOCATOR_H_

#include <stdexcept>
#include <algorithm>
#include <omp.h>



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
        char * memory_pool;
        size_t memory_pool_size;
        size_t volatile start_empty;
        size_t volatile start_full;
        size_t volatile full_blocks_vector_capacity;
        size_t volatile full_blocks_vector_size;
        mem_block /*volatile*/ * volatile full_blocks_vector_data;
        size_t transaction_end_fail;
        omp_lock_t lock_vector;
        omp_lock_t lock_transaction;
    public:
        cbmba_resource(void * memory_pool_, size_t memory_pool_size_)
            : memory_pool(reinterpret_cast<char*>(memory_pool_))
            , memory_pool_size(memory_pool_size_)
            , start_empty(0)
            , start_full(0)
            , full_blocks_vector_capacity(0)
            , full_blocks_vector_size(0)
            , full_blocks_vector_data(nullptr)
            , transaction_end_fail(~(size_t)0)
        {
            omp_init_lock(&lock_vector);
            omp_init_lock(&lock_transaction);
        }
        ~cbmba_resource()
        {
            omp_destroy_lock(&lock_vector);
            omp_destroy_lock(&lock_transaction);
        }
        cbmba_resource() = delete;
        cbmba_resource(cbmba_resource const &) = delete;
        cbmba_resource(cbmba_resource &&) = default;
        cbmba_resource & operator=(cbmba_resource const &) = delete;
        cbmba_resource & operator=(cbmba_resource &&) = default;
        void * allocate(size_t num_bytes, size_t align = 1)
        {
            if(num_bytes == 0) num_bytes = 1;
            if(num_bytes > memory_pool_size) throw std::runtime_error("Not enough memory in the pool");
            if(memory_pool_size % align != 0) throw std::runtime_error("Align has to divide memory_pool_size");
            if(reinterpret_cast<std::uintptr_t>(memory_pool) % align != 0) throw std::runtime_error("Memory pool has to be aligned the same way as is the requested align");

            size_t new_block_start;
            size_t new_block_start_aligned;
            size_t new_block_end;
            {
                omp_set_lock(&lock_vector);
                new_block_start = start_empty;
                new_block_start_aligned = ((new_block_start - 1) / align + 1) * align;
                new_block_end = new_block_start_aligned + num_bytes;
                if((new_block_end - 1) / memory_pool_size > new_block_start / memory_pool_size) // if the block would cross the pool size boundary
                {
                    new_block_start = ((start_empty - 1) / memory_pool_size + 1) * memory_pool_size;
                    new_block_start_aligned = ((new_block_start - 1) / align + 1) * align;
                    new_block_end = new_block_start_aligned + num_bytes;
                }

                add_to_full_blocks_vector(mem_block(new_block_start, new_block_start_aligned));
                start_empty = new_block_end;

                if(start_empty > transaction_end_fail) throw std::runtime_error("More memory requested in a single transaction then what is available in the pool");
                omp_unset_lock(&lock_vector);
            }

            // if(new_block_end > start_full + memory_pool_size) printf("Memory full, waiting. this=%p, ptr=%p\n", this, memory_pool + new_block_start_aligned % memory_pool_size);
            while(new_block_end > start_full + memory_pool_size) ; // busy waiting for some deallocation that moves the start_full pointer so I can return the memory to the user
            // printf("Memory granted, %zu B, this=%p, pool=%p--%p, ptr=%p\n", num_bytes, this, memory_pool, memory_pool + memory_pool_size, memory_pool + new_block_start_aligned % memory_pool_size);

            return memory_pool + new_block_start_aligned % memory_pool_size;
        }
        void deallocate(void * ptr)
        {
            if(ptr == nullptr) return;
            size_t mem_block_start_aligned = reinterpret_cast<char*>(ptr) - memory_pool;

            {
                omp_set_lock(&lock_vector);
                // printf("Deallocating, this=%p, ptr=%p\n", this, ptr);
                remove_from_full_blocks_vector(mem_block_start_aligned);
                if(full_blocks_vector_size == 0) start_full = start_empty = 0;
                else start_full = full_blocks_vector_data[0].start;
                omp_unset_lock(&lock_vector);
            }
        }
        void print_full_blocks()
        {
            printf("cbmba_resource full blocks:");
            for(size_t i = 0; i < full_blocks_vector_size; i++) printf(" %zu/%zu", full_blocks_vector_data[i].start, full_blocks_vector_data[i].start % memory_pool_size);
            printf("\n");
        }
        template<typename C> void do_transaction(const C & c)
        {
            // c must have an operator() with no parameters and must only do allocations using this resource
            omp_set_lock(&lock_transaction);
            transaction_end_fail = start_empty + memory_pool_size;
            c();
            transaction_end_fail = ~(size_t)0;
            omp_unset_lock(&lock_transaction);
        }
    private:
        void add_to_full_blocks_vector(mem_block mb)
        {
            if(full_blocks_vector_size == full_blocks_vector_capacity)
            {
                full_blocks_vector_capacity = 2 * full_blocks_vector_capacity + 1;
                mem_block * new_data = reinterpret_cast<mem_block*>(malloc(full_blocks_vector_capacity * sizeof(mem_block)));
                std::uninitialized_move_n(full_blocks_vector_data, full_blocks_vector_size, new_data);
                free(full_blocks_vector_data);
                full_blocks_vector_data = new_data;
            }
            std::uninitialized_move_n(&mb, 1, full_blocks_vector_data + full_blocks_vector_size);
            full_blocks_vector_size++;
        }
        void remove_from_full_blocks_vector(size_t mem_block_start_aligned)
        {
            mem_block * end_iter = full_blocks_vector_data + full_blocks_vector_size;
            mem_block * it = std::find_if(full_blocks_vector_data, end_iter, [=](mem_block const & mb){return mb.start_aligned % memory_pool_size == mem_block_start_aligned;});
            if(it == end_iter) throw std::runtime_error("Bad remove from full blocks vector");
            std::move(it + 1, end_iter, it);
            full_blocks_vector_size--;
        }
    };



    template<bool ad, bool ah>
    class cbmba
    {
    public:
        static constexpr bool is_data_device_accessible = ad;
        static constexpr bool is_data_host_accessible = ah;
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
