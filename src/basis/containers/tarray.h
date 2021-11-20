
#ifndef SRC_BASIS_CONTAINERS_TARRAY_H_
#define SRC_BASIS_CONTAINERS_TARRAY_H_

#include <cmath>
#include <vector>
#include <cstring>

#include <omp.h>

namespace espreso {


/// Threaded array
template <typename TType>
class tarray {

public:
	static std::vector<TType> distribute(int threads, TType size, TType minchunk = 0);

	tarray(size_t threads, size_t size, TType init = TType{});
	tarray(size_t threads, const std::vector<TType> &data);
	tarray(const std::vector<std::vector<TType> > &data);
	tarray(const std::vector<TType> &data);
	tarray(const std::vector<size_t> &distribution, const std::vector<TType> &data);
	template <typename TOther>
	tarray(const std::vector<size_t> &distribution, TOther begin, TOther end);
	tarray(const std::vector<size_t> &distribution, size_t duplication, TType init = TType{});

	tarray(const tarray<TType> &other);
	tarray(tarray<TType> &&other);
	tarray<TType>& operator=(const tarray<TType> &other);
	tarray<TType>& operator=(tarray<TType> &&other);

	size_t                     size()         const { return _size; }
	size_t                     size(size_t t) const { return _distribution[t + 1] - _distribution[t]; }
	size_t                     threads()      const { return _distribution.size() - 1; }
	const std::vector<size_t>& distribution() const { return _distribution; }

	TType& operator[] (size_t n) { return _data[n]; }
	TType* data()                { return _data; }
	TType& back()                { return _data[_size - 1]; }
	TType& front()               { return _data[0]; }
	TType* begin()               { return _data; }
	TType* end()                 { return _data + _size; }
	TType* begin(size_t t)       { return _data + _distribution[t]; }
	TType* end(size_t t)         { return _data + _distribution[t + 1]; }

	const TType& operator[] (size_t n) const { return _data[n]; }
	const TType* data()                const { return _data; }
	const TType& back()                const { return _data[_size - 1]; }
	const TType& front()               const { return _data[0]; }
	const TType* begin()               const { return _data; }
	const TType* end()                 const { return _data + _size; }
	const TType* cbegin()              const { return _data; }
	const TType* cend()                const { return _data + _size; }
	const TType* cbegin(size_t t)      const { return _data + _distribution[t]; }
	const TType* cend(size_t t)        const { return _data + _distribution[t + 1]; }

	size_t packedSize() const;
	void pack(char* &p) const;
	void unpack(const char* &p);

	~tarray();

private:
	size_t _size;
	TType *_data;

	std::vector<size_t> _distribution;
};

}

#include "tarray.hpp"

#endif /* SRC_BASIS_CONTAINERS_TARRAY_H_ */
