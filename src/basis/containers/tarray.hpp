
#ifndef SRC_BASIS_CONTAINERS_TARRAY_HPP_
#define SRC_BASIS_CONTAINERS_TARRAY_HPP_

#include "tarray.h"

namespace espreso {

template <typename TType>
std::vector<TType> tarray<TType>::distribute(int threads, TType size, TType minchunk)
{
	std::vector<TType> distribution(threads + 1);

	TType chunksize = std::max((TType)std::ceil(size / (double)threads), minchunk);
	for (int t = 1; t < threads; t++) {
		distribution[t] = t * chunksize;
		if (distribution[t] > size) {
			distribution[t] = size;
		}
	}
	distribution[threads] = size;

	return distribution;
}

template <typename TType>
tarray<TType>::tarray(const tarray<TType> &other)
: _size(other._size), _data(other._data), _distribution(other._distribution)
{
	if (_size) {
		_data = new TType[_size];
		#pragma omp parallel for
		for (size_t t = 0; t < threads(); t++) {
			memcpy(reinterpret_cast<void*>(_data + _distribution[t]), other._data + _distribution[t], (_distribution[t + 1] - _distribution[t]) * sizeof(TType));
		}
	}
}

template <typename TType>
tarray<TType>::tarray(tarray<TType> &&other)
: _size(std::move(other._size)), _data(std::move(other._data)), _distribution(std::move(other._distribution))
{
	other._size = 0;
	other._data = NULL;
}

template <typename TType>
tarray<TType>& tarray<TType>::operator=(const tarray<TType> &other)
{
	if (this != &other) {
		_size = other._size;
		_distribution = other._distribution;

		if (_data) { delete[] _data; }
		_data = new TType[_size];
		#pragma omp parallel for
		for (size_t t = 0; t < other.threads(); t++) {
			memcpy(_data + _distribution[t], other._data + _distribution[t], (_distribution[t + 1] - _distribution[t]) * sizeof(TType));
		}
	}
	return *this;
}

template <typename TType>
tarray<TType>& tarray<TType>::operator=(tarray<TType> &&other)
{
	if (this != &other) {
		_size = std::move(other._size);
		if (_data) { delete[] _data; }
		_data = std::move(other._data);
		_distribution = std::move(other._distribution);

		other._size = 0;
		other._data = NULL;
	}
	return *this;
}

template <typename TType>
tarray<TType>::tarray(size_t threads, size_t size, TType init)
:  _size(size), _data(NULL)
{
	_distribution = tarray<size_t>::distribute(threads, size);

	if (_size) {
		_data = new TType[_size];
		#pragma omp parallel for
		for (size_t t = 0; t < threads; t++) {
			for (size_t i = _distribution[t]; i < _distribution[t + 1]; i++) {
				_data[i] = init;
			}
		}
	}
}

template <typename TType>
tarray<TType>::tarray(size_t threads, const std::vector<TType> &data)
: _size(data.size()), _data(NULL), _distribution(tarray<size_t>::distribute(threads, data.size()))
{
	if (_size) {
		_data = new TType[_size];
		#pragma omp parallel for
		for (size_t t = 1; t < _distribution.size(); t++) {
			memcpy(reinterpret_cast<void*>(_data + _distribution[t - 1]), data.data() + _distribution[t - 1], (_distribution[t] - _distribution[t - 1]) * sizeof(TType));
		}
	}
}

template <typename TType>
tarray<TType>::tarray(const std::vector<std::vector<TType> > &data)
: _size(0), _data(NULL)
{
	_distribution = std::vector<size_t>(data.size() + 1, 0);
	for (size_t t = 0; t < data.size(); t++) {
		_size += data[t].size();
		_distribution[t + 1] = _size;
	}

	if (_size) {
		_data = new TType[_size];
		#pragma omp parallel for
		for (size_t t = 0; t < threads(); t++) {
			memcpy(reinterpret_cast<void*>(_data + _distribution[t]), data[t].data(), data[t].size() * sizeof(TType));
		}
	}
}

template <typename TType>
tarray<TType>::tarray(const std::vector<TType> &data)
: _data(NULL)
{
	_size = data.size();
	_distribution = { 0, _size };

	if (_size) {
		_data = new TType[_size];
		memcpy(reinterpret_cast<void*>(_data), data.data(), data.size() * sizeof(TType));
	}
}

template <typename TType>
tarray<TType>::tarray(const std::vector<size_t> &distribution, const std::vector<TType> &data)
: _size(data.size()), _data(NULL), _distribution(distribution)
{
	if (_size) {
		_data = new TType[_size];
		#pragma omp parallel for
		for (size_t t = 1; t < _distribution.size(); t++) {
			memcpy(reinterpret_cast<void*>(_data + _distribution[t - 1]), data.data() + _distribution[t - 1], (_distribution[t] - _distribution[t - 1]) * sizeof(TType));
		}
	}
}

template <typename TType>
template <typename TOther>
tarray<TType>::tarray(const std::vector<size_t> &distribution, TOther begin, TOther end)
: _size(end - begin), _data(NULL), _distribution(distribution)
{
	if (_size) {
		_data = new TType[_size];
		#pragma omp parallel for
		for (size_t t = 1; t < _distribution.size(); t++) {
			std::copy(begin + _distribution[t - 1], begin + _distribution[t], _data + _distribution[t - 1]);
		}
	}
}

template <typename TType>
tarray<TType>::tarray(const std::vector<size_t> &distribution, size_t duplication, TType init)
: _size(duplication * distribution.back()), _data(NULL), _distribution(distribution)
{
	for (size_t t = 1; t < _distribution.size(); t++) {
		_distribution[t] *= duplication;
	}

	if (_size) {
		_data = new TType[_size];
		#pragma omp parallel for
		for (size_t t = 0; t < _distribution.size() - 1; t++) {
			for (size_t i = _distribution[t]; i < _distribution[t + 1]; i++) {
				_data[i] = init;
			}
		}
	}
}

template <typename TType>
size_t tarray<TType>::packedSize() const
{
	return
			sizeof(size_t) + _distribution.size() * sizeof(size_t) +
			sizeof(size_t) +
			_size * sizeof(TType);
}

template <typename TType>
void tarray<TType>::pack(char* &p) const
{
	size_t dsize = _distribution.size();
	memcpy(p, &dsize, sizeof(size_t));
	p += sizeof(size_t);
	memcpy(p, _distribution.data(), _distribution.size() * sizeof(size_t));
	p += _distribution.size() * sizeof(size_t);
	memcpy(p, &_size, sizeof(size_t));
	p += sizeof(size_t);
	if (_size) {
		memcpy(p, _data, _size * sizeof(TType));
		p += _size * sizeof(TType);
	}
}

template <typename TType>
void tarray<TType>::unpack(const char* &p)
{
	size_t dsize;
	memcpy(&dsize, p, sizeof(size_t));
	p += sizeof(size_t);
	_distribution.resize(dsize);
	memcpy(_distribution.data(), p, _distribution.size() * sizeof(size_t));
	p += _distribution.size() * sizeof(size_t);

	memcpy(&_size, p, sizeof(size_t));
	p += sizeof(size_t);

	if (_data != NULL && _distribution.back() != _size) {
		delete _data;
		_data = NULL;
	}
	if (_size) {
		if (_data == NULL) {
			_data = new TType[_size];
		}
		memcpy(reinterpret_cast<void*>(_data), p, _size * sizeof(TType));
		p += _size * sizeof(TType);
	}
}

template <typename TType>
tarray<TType>::~tarray()
{
	if (_data) {
		delete[] _data;
	}
}

}

#endif /* SRC_BASIS_CONTAINERS_TARRAY_HPP_ */
