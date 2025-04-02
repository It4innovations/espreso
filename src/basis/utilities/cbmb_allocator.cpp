
#include "cbmb_allocator.h"

#include "esinfo/eslog.hpp"

#include <omp.h>



namespace espreso {

    cbmba_resource::cbmba_resource(void * memory_pool_, size_t memory_pool_size_)
        : memory_pool(reinterpret_cast<char*>(memory_pool_))
        , memory_pool_size(memory_pool_size_)
        , start_empty(0)
        , start_full(0)
        , transaction_end_fail(~(size_t)0)
    {
    }

    cbmba_resource::~cbmba_resource()
    {
    }

    void * cbmba_resource::allocate(size_t num_bytes, size_t align)
    {
        if(num_bytes == 0) return nullptr;
        if(num_bytes > memory_pool_size) eslog::error("CBMBA: single allocation cannot be larger than capacity\n");
        if(memory_pool_size % align != 0) eslog::error("CBMBA: align has to divide memory_pool_size\n");
        if(reinterpret_cast<std::uintptr_t>(memory_pool) % align != 0) eslog::error("CBMBA: memory pool has to be aligned the same way as is the requested align\n");

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

            if(in_transaction && start_empty > transaction_end_fail) eslog::error("CBMBA: more memory requested in a single transaction then what is available in the pool\n");

            if(in_transaction) mem_allocated_in_transaction += num_bytes;
            if(in_transaction && mem_allocated_in_transaction > memory_pool_size / 2) do_reset_to_start_after_transaction = true;

            // eslog::info("Waiting for memory, %zu B, this=%p, new_block_end=%zu, sf+mps=%zu\n", num_bytes, this, new_block_end, start_full + memory_pool_size);
            cv.wait(lk, [&]{ return new_block_end <= start_full + memory_pool_size; });
            // eslog::info("Memory granted,     %zu B, this=%p, new_block_end=%zu, sf+mps=%zu\n", num_bytes, this, new_block_end, start_full + memory_pool_size);
        }

        return memory_pool + new_block_start_aligned % memory_pool_size;
    }

    void cbmba_resource::deallocate(void * & ptr)
    {
        if(ptr == nullptr) return;
        size_t mem_block_start_aligned = reinterpret_cast<char*>(ptr) - memory_pool;

        {
            std::lock_guard<std::mutex> lk(mtx_blocks);

            // eslog::info("Deallocating, this=%p, ptr=%p, pos=%zu\n", this, ptr, mem_block_start_aligned);
            // remove_from_full_blocks_vector(mem_block_start_aligned);
            auto it = std::find_if(full_blocks.begin(), full_blocks.end(), [=](mem_block const & mb){return mb.start_aligned % memory_pool_size == mem_block_start_aligned;});
            if(it == full_blocks.end()) eslog::error("CBMBA: bad remove from full blocks vector\n");
            full_blocks.erase(it);

            if(full_blocks.size() == 0) start_full = start_empty = 0;
            else start_full = full_blocks.front().start;
        }
        cv.notify_all();
        ptr = nullptr;
    }

    void cbmba_resource::print_full_blocks()
    {
        eslog::info("cbmba_resource full blocks:");
        for(mem_block & mb : full_blocks) eslog::info(" %zu/%zu", mb.start, mb.start % memory_pool_size);
        eslog::info("\n");
    }

    void cbmba_resource::do_transaction(const std::function<void(void)> & f)
    {
        std::lock_guard<std::mutex> lk(mtx_transaction);
        // f must do allocations using only this resource
        transaction_end_fail = start_empty + memory_pool_size;
        mem_allocated_in_transaction = 0;
        in_transaction = true;
        f();
        in_transaction = false;
        if(do_reset_to_start_after_transaction) start_empty = ((start_empty - 1) / memory_pool_size + 1) * memory_pool_size;
        mem_allocated_in_transaction = 0;
    }

    size_t cbmba_resource::get_curr_memory_requested_in_transaction() const
    {
        return mem_allocated_in_transaction;
    }

    size_t cbmba_resource::get_max_capacity() const
    {
        return memory_pool_size;
    }

}
