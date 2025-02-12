
#ifndef SRC_BASIS_CONTAINERS_ALLOCATORS_H_
#define SRC_BASIS_CONTAINERS_ALLOCATORS_H_

#include <memory>
#include <cstdlib>

namespace espreso {

// allocator that allows to mimic new[] behavior
// it skipes initialization when vector.resize(size) is called
template<typename T, typename A = std::allocator<T> >
class initless_allocator: public A {
    typedef std::allocator_traits<A> a_t;
public:
    template<typename U> struct rebind {
        using other = initless_allocator<U, typename a_t::template rebind_alloc<U> >;
    };

    using A::A;

    template<typename U>
    void construct(U* ptr) noexcept(std::is_nothrow_default_constructible<U>::value) {
        ::new (static_cast<void*>(ptr)) U;
    }
    template<typename U, typename ...Args>
    void construct(U* ptr, Args&&... args) {
        a_t::construct(static_cast<A&>(*this), ptr, std::forward<Args>(args)...);
    }
};

// not compatible with std::allocator
class cpu_allocator
{
public:
    static constexpr bool is_data_host_accessible = true;
    static constexpr bool is_data_device_accessible = false;
    static constexpr bool always_equal = true;
public:
    void * allocate(size_t num_bytes, size_t alignment)
    {
        return aligned_alloc(alignment, num_bytes);
    }
    void * allocate(size_t num_bytes)
    {
        return allocate(num_bytes, 1);
    }
    template<typename T>
    T * allocate(size_t count)
    {
        return reinterpret_cast<T*>(allocate(count * sizeof(T), sizeof(T)));
    }
    template<typename T>
    void deallocate(T * ptr)
    { 
        free(ptr);
    }
    bool operator==(const cpu_allocator & /*other*/)
    {
        return true;
    }
    bool operator!=(const cpu_allocator & other)
    {
        return !(*this == other);
    }
};

}

#endif /* SRC_BASIS_CONTAINERS_ALLOCATORS_H_ */
