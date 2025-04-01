
#ifndef SRC_BASIS_UTILITIES_ARENA_ALLOCATOR_H_
#define SRC_BASIS_UTILITIES_ARENA_ALLOCATOR_H_

#include <mutex>
#include <memory>



namespace espreso {

    class arena_allocator_resource
    {
    private:
        char * mem;
        size_t capacity;
        size_t curr_offset = 0;
        size_t align_B;
        std::mutex mtx;
        size_t num_allocations = 0;
    public:
        arena_allocator_resource() = delete;
        arena_allocator_resource(arena_allocator_resource const &) = delete;
        arena_allocator_resource(arena_allocator_resource &&) = delete;
        ~arena_allocator_resource();
        arena_allocator_resource & operator=(arena_allocator_resource const &) = delete;
        arena_allocator_resource & operator=(arena_allocator_resource &&) = delete;
    public:
        arena_allocator_resource(void * mem_, size_t capacity_, size_t align_B_);
        void * allocate(size_t num_bytes);
        template<typename T>
        T * allocate(size_t count)
        {
            return reinterpret_cast<T*>(allocate(count * sizeof(T)));
        }
        void deallocate(void * & ptr);
        template<typename T>
        void deallocate(T * & ptr)
        {
            deallocate((void*&)ptr);
        }
        void reset_arena();
        size_t get_free_size();
        size_t get_capacity();
    };

    template<bool ad, bool ah>
    class arena_allocator
    {
    public:
        static constexpr bool is_data_device_accessible = ad;
        static constexpr bool is_data_host_accessible = ah;
        static constexpr bool always_equal = false;
    private:
        std::shared_ptr<arena_allocator_resource> internal;
    public:
        arena_allocator() = delete;
        arena_allocator(arena_allocator<ad,ah> const &) = default;
        arena_allocator(arena_allocator<ad,ah> &&) = default;
        ~arena_allocator() = default;
        arena_allocator<ad,ah> & operator=(arena_allocator<ad,ah> const &) = default;
        arena_allocator<ad,ah> & operator=(arena_allocator<ad,ah> &&) = default;
    public:
        arena_allocator(void * mem, size_t capacity, size_t align_B)
        {
            internal = std::make_unique<arena_allocator_resource>(mem, capacity, align_B);
        }
        void * allocate(size_t num_bytes)
        {
            return internal->allocate(num_bytes);
        }
        template<typename T>
        T * allocate(size_t count)
        {
            return reinterpret_cast<T*>(allocate(count * sizeof(T)));
        }
        template<typename T>
        void deallocate(T * & ptr)
        {
            internal->deallocate(ptr);
        }
        void reset_arena()
        {
            internal->reset_arena();
        }
        size_t get_free_size()
        {
            return internal->get_free_size();
        }
        size_t get_capacity()
        {
            return internal->get_capacity();
        }
    };

    using arena_d = arena_allocator<true,false>;
    using arena_h = arena_allocator<false,true>;

}

#endif /* SRC_BASIS_UTILITIES_ARENA_ALLOCATOR_H_ */
