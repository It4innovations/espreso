
#ifndef SRC_BASIS_CONTAINERS_ALLOCATORS_H_
#define SRC_BASIS_CONTAINERS_ALLOCATORS_H_

#include <memory>

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

}

#endif /* SRC_BASIS_CONTAINERS_ALLOCATORS_H_ */
