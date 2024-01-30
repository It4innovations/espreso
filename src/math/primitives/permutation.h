
#ifndef SRC_MATH2_PRIMITIVES_PERMUTATION_H_
#define SRC_MATH2_PRIMITIVES_PERMUTATION_H_

#include <algorithm>

#include "basis/containers/allocators.h"
#include "esinfo/eslog.hpp"

namespace espreso {

template <typename I>
class _Permutation {
	public:
    I * dst_to_src = nullptr;   // idx_src = dst_to_src[idx_dst];   dst[i] = src[dst_to_src[i]];
    I * src_to_dst = nullptr;   // idx_dst = src_to_dst[idx_src];   dst[src_to_dst[i]] = src[i];
	I size;
};

template <typename I = int, typename A = cpu_allocator>
class Permutation: public _Permutation<I>
{
public:
	Permutation(const A &ator_ = A()): _Permutation<I>{}, _allocated{}, ator(ator_)
	{

	}

	Permutation(const Permutation &other): _Permutation<I>{}, _allocated{}, ator(other.ator)
	{
		static_assert(A::is_data_host_accessible, "the allocator does not provide host-accessible memory\n");
		realloc(_allocated, other.size);
		_Permutation<I>::operator=(_allocated);
		std::copy_n(other.dst_to_src, other.size, this->dst_to_src);
		std::copy_n(other.src_to_dst, other.size, this->src_to_dst);
	}

	Permutation(Permutation &&other): _Permutation<I>{}, _allocated{}, ator(std::move(other.ator))
	{
		swap(*this, other);
		swap(_allocated, other._allocated);
	}

	Permutation& operator=(const Permutation &other) = delete;

	Permutation& operator=(Permutation &&other) = delete;

	~Permutation()
	{
		clear(_allocated);
	}

	void resize(I size)
	{
		realloc(_allocated, size);
		_Permutation<I>::operator=(_allocated);
	}

	template<typename I2, typename A2>
	void resize(const Permutation<I2,A2> &other)
	{
		resize(other.size);
	}

	void pattern(const Permutation &other)
	{
		if constexpr(!A::always_equal) if(this->ator != other.ator) eslog::error("not implemented for unequal allocators\n");
		realloc(_allocated, other.size);
		_Permutation<I>::operator=(_allocated);
	}

	void shallowCopy(const Permutation &other)
	{
		if constexpr(!A::always_equal) if(this->ator != other.ator) eslog::error("not implemented for unequal allocators\n");
		_Permutation<I>::operator=(other);
	}

	void swap(Permutation &other)
	{
		swap(*this, other);
		swap(_allocated, other._allocated);
		std::swap(this->ator, other.ator);
	}

	void invert(I * map_in, I * map_out)
	{
		static_assert(A::is_data_host_accessible, "the allocator does not provide host-accessible memory\n");
		for(I i = 0; i < this->size; i++)
		{
			map_out[map_in[i]] = i;
		}
	}

	void clear()
	{
		clear(_allocated);
		resize(0);
	}

	_Permutation<I> _allocated;
	A ator;

protected:
	template <typename Type>
	void _swap(Type &p, Type &q)
	{
		Type tmp = p; p = q; q = tmp;
	}

	void swap(_Permutation<I> &p, _Permutation<I> &q)
	{
		_swap(p.size, q.size);
		_swap(p.vals, q.vals);
	}

	void realloc(_Permutation<I> &p, I newsize)
	{
		if (p.size < newsize) {
			clear(p);
			p.dst_to_src = ator.template allocate<I>(newsize);
			p.src_to_dst = ator.template allocate<I>(newsize);
		}
		p.size = newsize;
	}

	void clear(_Permutation<I> &p)
	{
		if (p.dst_to_src) { ator.deallocate(p.dst_to_src); p.dst_to_src = nullptr; }
		if (p.src_to_dst) { ator.deallocate(p.src_to_dst); p.src_to_dst = nullptr; }
	}
};

}

#endif /* SRC_MATH2_PRIMITIVES_PERMUTATION_H_ */
