
#ifndef SRC_BASIS_CONTAINERS_TARRAY_HPP_
#define SRC_BASIS_CONTAINERS_TARRAY_HPP_

#include "tarray.h"
#include <memory>

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
: _size(other._size), _alignment(other._alignment), _data(NULL), _dataUnaligned(NULL), _distribution(other._distribution)
{
	if (_size) {
		this->allocate();
		#pragma omp parallel for
		for (size_t t = 0; t < threads(); t++) {
			memcpy(reinterpret_cast<void*>(_data + _distribution[t]), other._data + _distribution[t], (_distribution[t + 1] - _distribution[t]) * sizeof(TType));
		}
	}
}

template <typename TType>
tarray<TType>::tarray(tarray<TType> &&other)
: _size(std::move(other._size)), _alignment(std::move(other._alignment)), _data(std::move(other._data)), _dataUnaligned(std::move(other._dataUnaligned)), _distribution(std::move(other._distribution))
{
	other._size = 0;
	other._dataUnaligned = NULL;
}

template <typename TType>
tarray<TType>& tarray<TType>::operator=(const tarray<TType> &other)
{
	if (this != &other) {
		_size = other._size;
		_distribution = other._distribution;
		_alignment = other._alignment;

		this->deallocate();
		this->allocate();
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
		this->deallocate();
		_dataUnaligned = std::move(other._dataUnaligned);
		_data = std::move(other._data);
		_distribution = std::move(other._distribution);
		_alignment = std::move(other._alignment);

		other._size = 0;
		other._dataUnaligned = NULL;
		other._data = NULL;
	}
	return *this;
}

template <typename TType>
tarray<TType>::tarray(size_t threads, size_t size, TType init, size_t alignment)
: _size(size), _alignment(alignment), _data(NULL), _dataUnaligned(NULL)
{
	_distribution = tarray<size_t>::distribute(threads, size);

	if (_size) {
		this->allocate();
		#pragma omp parallel for
		for (size_t t = 0; t < threads; t++) {
			for (size_t i = _distribution[t]; i < _distribution[t + 1]; i++) {
				_data[i] = init;
			}
		}
	}
}

template <typename TType>
tarray<TType>::tarray(size_t threads, size_t size, bool skipinit, size_t alignment)
: _size(size), _alignment(alignment), _data(NULL), _dataUnaligned(NULL)
{
	_distribution = tarray<size_t>::distribute(threads, size);

	if (_size) {
		this->allocate();
		if (!skipinit) {
			#pragma omp parallel for
			for (size_t t = 0; t < threads; t++) {
				for (size_t i = _distribution[t]; i < _distribution[t + 1]; i++) {
					_data[i] = TType{};
				}
			}
		}
	}
}

template <typename TType>
tarray<TType>::tarray(size_t threads, const std::vector<TType> &data, size_t alignment)
: _size(data.size()), _alignment(alignment), _data(NULL),  _dataUnaligned(NULL), _distribution(tarray<size_t>::distribute(threads, data.size()))
{
	if (_size) {
		this->allocate();
		#pragma omp parallel for
		for (size_t t = 1; t < _distribution.size(); t++) {
			memcpy(reinterpret_cast<void*>(_data + _distribution[t - 1]), data.data() + _distribution[t - 1], (_distribution[t] - _distribution[t - 1]) * sizeof(TType));
		}
	}
}

template <typename TType>
tarray<TType>::tarray(const std::vector<std::vector<TType> > &data, size_t alignment)
: _size(0), _alignment(alignment), _data(NULL),  _dataUnaligned(NULL)
{
	_distribution = std::vector<size_t>(data.size() + 1, 0);
	for (size_t t = 0; t < data.size(); t++) {
		_size += data[t].size();
		_distribution[t + 1] = _size;
	}

	if (_size) {
		this->allocate();
		#pragma omp parallel for
		for (size_t t = 0; t < threads(); t++) {
			memcpy(reinterpret_cast<void*>(_data + _distribution[t]), data[t].data(), data[t].size() * sizeof(TType));
		}
	}
}

template <typename TType>
tarray<TType>::tarray(const std::vector<TType> &data, size_t alignment)
: _alignment(alignment), _data(NULL), _dataUnaligned(NULL)
{
	_size = data.size();
	_distribution = { 0, _size };

	if (_size) {
		this->allocate();
		memcpy(reinterpret_cast<void*>(_data), data.data(), data.size() * sizeof(TType));
	}
}

template <typename TType>
tarray<TType>::tarray(const std::vector<size_t> &distribution, const std::vector<TType> &data, size_t alignment)
: _size(data.size()), _alignment(alignment), _dataUnaligned(NULL), _distribution(distribution)
{
	if (_size) {
		this->allocate();
		#pragma omp parallel for
		for (size_t t = 1; t < _distribution.size(); t++) {
			memcpy(reinterpret_cast<void*>(_data + _distribution[t - 1]), data.data() + _distribution[t - 1], (_distribution[t] - _distribution[t - 1]) * sizeof(TType));
		}
	}
}

template <typename TType>
tarray<TType>::tarray(const std::vector<size_t> &distribution, size_t duplication, TType init, size_t alignment)
: _size(duplication * distribution.back()), _alignment(alignment), _data(NULL), _dataUnaligned(NULL), _distribution(distribution)
{
	for (size_t t = 1; t < _distribution.size(); t++) {
		_distribution[t] *= duplication;
	}

	if (_size) {
		this->allocate();
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
		memcpy(p, _dataUnaligned, _size * sizeof(TType));
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

	if (_dataUnaligned != NULL && _distribution.back() != _size) {
		this->deallocate();
	}
	if (_size) {
		this->allocate();
		memcpy(reinterpret_cast<void*>(_dataUnaligned), p, _size * sizeof(TType));
		p += _size * sizeof(TType);
	}
}

template <typename TType>
tarray<TType>::~tarray()
{
	this->deallocate();
}

template <typename TType>
void tarray<TType>::allocate()
{
	if(_alignment != 0){
		size_t sizeInElements = _size + ((_alignment - 1) / sizeof(TType)) +1;
		size_t sizeInBytes = sizeInElements*sizeof(TType);
		_dataUnaligned = new TType[sizeInElements];
		void* tmp = static_cast<void *>(_dataUnaligned);
		_data = static_cast<TType *>(std::align(_alignment, _size*sizeof(TType), tmp, sizeInBytes));
	}
	else
	{
		_dataUnaligned = new TType[_size];
		_data = _dataUnaligned;
	}
}

template <typename TType>
void tarray<TType>::deallocate()
{
	if (_dataUnaligned) {
		delete[] _dataUnaligned;
		_dataUnaligned = NULL;
		_data = NULL;
	}
}

}

#endif /* SRC_BASIS_CONTAINERS_TARRAY_HPP_ */
