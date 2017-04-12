
#ifndef SRC_BASIS_CONTAINERS_PARRAY_H_
#define SRC_BASIS_CONTAINERS_PARRAY_H_

#include <cstddef>

namespace espreso {

template <typename TType>
class parray {

public:
	static size_t* ditribute(size_t threads, size_t size);

	parray(const parray<TType> &other);
	parray(parray<TType> &&other);
	parray(size_t threads, size_t size);
	parray(size_t threads, const TType *data, size_t size);
	parray(size_t threads, const std::vector<TType> &data);
	parray(const std::vector<std::vector<TType> > &data);

	parray<TType>& operator=(const parray<TType> &other);

	const size_t* distribution() const { return _distribution; }
	size_t threads() const { return _threads; }
	size_t size() const { return _size; }

	TType& operator[] (size_t n) { return _data[n]; }
	TType* data() { return _data; }
	TType& back() { return _data[_size - 1]; }
	TType& front() { return _data[0]; }
	TType* begin() { return _data; }
	TType* end() { return _data + _size; }

	const TType& operator[] (size_t n) const { return _data[n]; }
	const TType* data() const { return _data; }
	const TType& back() const { return _data[_size - 1]; }
	const TType& front() const { return _data[0]; }
	const TType* begin() const { return _data; }
	const TType* end() const { return _data + _size; }
	const TType* cbegin() const { return _data; }
	const TType* cend() const { return _data + _size; }

	~parray();

private:
	size_t _threads;
	size_t _size;
	size_t *_distribution;
	TType *_data;
};

}

#include "../containers/parray.hpp"



#endif /* SRC_BASIS_CONTAINERS_PARRAY_H_ */
