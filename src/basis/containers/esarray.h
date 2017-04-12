
#ifndef SRC_BASIS_CONTAINERS_ESARRAY_H_
#define SRC_BASIS_CONTAINERS_ESARRAY_H_

#include <cstddef>

namespace espreso {

/**
 * Simple structure of data pointer and its size.
 */

template <typename TType>
class esarray {

public:
	esarray(const esarray<TType> &other): _size(other._size), _data(other._data) {}
	esarray(esarray<TType> &&other): _size(std::move(other._size)), _data(std::move(other._data)) {}
	esarray(size_t size, TType *data): _size(size), _data(data) {}

	void reinit(size_t size, TType *data) { _size = size; _data = data; }

	esarray<TType>& operator=(const esarray<TType> &other)
	{
		if (this != &other) {
			_size = other._size;
			_data = other._data;
		}
		return *this;
	}

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

private:
	size_t _size;
	TType *_data;
};

}



#endif /* SRC_BASIS_CONTAINERS_ESARRAY_H_ */
