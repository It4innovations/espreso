
#ifndef SRC_BASIS_CONTAINERS_PARRAY_HPP_
#define SRC_BASIS_CONTAINERS_PARRAY_HPP_

#include "../containers/parray.h"

#include <cmath>

namespace espreso {

template <typename TType>
size_t* parray<TType>::ditribute(size_t threads, size_t size)
{
	size_t *distribution = new size_t[threads + 1];

	size_t chunkSize = std::ceil(size / (double)threads);
	for (size_t t = 1; t < threads; t++) {
		distribution[t] = t * chunkSize;
		if (distribution[t] > size) {
			distribution[t] = size;
		}
	}
	distribution[threads] = size;

	return distribution;
}


template <typename TType>
parray<TType>::parray(const parray<TType> &other)
: _threads(other._threads), _size(other._size), _distribution(other._distribution)
{
	_data = new TType[_size];
	#pragma omp parallel for
	for (size_t t = 0; t < threads; t++) {
		for (size_t i = _distribution[t]; i < _distribution[t + 1]; i++) {
			_data[i] = other._data[i];
		}
	}
}

template <typename TType>
parray<TType>::parray(parray<TType> &&other)
: _threads(std::move(other._threads)), _size(std::move(other._size)), _distribution(std::move(other._distribution)), _data(std::move(other._data))
{

}

template <typename TType>
parray<TType>::parray(size_t threads, size_t size)
: _threads(threads), _size(size), _distribution(ditribute(threads, size))
{
	_data = new TType[_size];
	#pragma omp parallel for
	for (size_t t = 0; t < threads; t++) {
		for (size_t i = _distribution[t]; i < _distribution[t + 1]; i++) {
			_data[i] = TType{};
		}
	}
}

template <typename TType>
parray<TType>::parray(size_t threads, const TType *data, size_t size)
: _threads(threads), _size(size), _distribution(ditribute(threads, size))
{
	_data = new TType[_size];
	#pragma omp parallel for
	for (size_t t = 0; t < threads; t++) {
		memcpy(_data + _distribution[t], data + _distribution[t], (_distribution[t + 1] - _distribution[t]) * sizeof(TType));
	}
}

template <typename TType>
parray<TType>::parray(size_t threads, const std::vector<TType> &data)
: _threads(threads), _size(data.size()), _distribution(ditribute(threads, data.size()))
{
	_data = new TType[_size];
	#pragma omp parallel for
	for (size_t t = 0; t < _threads; t++) {
		memcpy(_data + _distribution[t], data.data() + _distribution[t], (_distribution[t + 1] - _distribution[t]) * sizeof(TType));
	}
}

template <typename TType>
parray<TType>::parray(const std::vector<std::vector<TType> > &data)
: _threads(data.size()), _size(0)
{
	_distribution = new size_t[data.size() + 1];
	_distribution[0] = 0;
	for (size_t t = 0; t < data.size(); t++) {
		_size += data[t].size();
		_distribution[t + 1] = _size;
	}

	_data = new TType[_size];
	#pragma omp parallel for
	for (size_t t = 0; t < _threads; t++) {
		memcpy(_data + _distribution[t], data[t].data(), data[t].size() * sizeof(TType));
	}
}

template <typename TType>
parray<TType>& parray<TType>::operator=(const parray<TType> &other)
{
	if (this != &other) {
		_threads = other._threads;
		_size = other._size;

		delete[] _distribution;
		_distribution = new size_t[_threads + 1];
		memcpy(_distribution, other._distribution, (_threads + 1) * sizeof(size_t));

		delete[] _data;
		_data = new TType[_size];
		#pragma omp parallel for
		for (size_t t = 0; t < _threads; t++) {
			memcpy(_data + _distribution[t], other._data + _distribution[t], (_distribution[t + 1] - _distribution[t]) * sizeof(TType));
		}
	}
	return *this;
}

template <typename TType>
parray<TType>::~parray()
{
	delete[] _distribution;
	delete[] _data;
}

}



#endif /* SRC_BASIS_CONTAINERS_PARRAY_HPP_ */
