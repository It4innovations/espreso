
#ifndef SRC_BASIS_CONTAINERS_PCRSDATA_HPP_
#define SRC_BASIS_CONTAINERS_PCRSDATA_HPP_

#include "pcrsdata.h"

namespace espreso {


template <typename TType>
pcrsdata<TType>::pcrsdata(parray<TType> &&crs, parray<TType> &&data)
: _crs(std::move(crs)), _data(std::move(data))
{
	_iterator = new crsiterator<TType>(0, _crs.size() - 1, _crs.data(), _data.data());
	_titerators = new crsiterator<TType>*[threads()];
	#pragma omp parallel for
	for (size_t t = 0; t < _data.threads() - 1; t++) {
		_titerators[t] = new crsiterator<TType>(_crs.distribution()[t], _crs.distribution()[t + 1], _crs.data(), _data.data());
	}
	size_t t = _data.threads() - 1;
	_titerators[t] = new crsiterator<TType>(_crs.distribution()[t], _crs.distribution()[t + 1] - 1, _crs.data(), _data.data());
}

template <typename TType>
static std::vector<std::vector<TType> >& push_to_last(std::vector<std::vector<TType> > &vector)
{
	vector.back().push_back(TType{});
	return vector;
}

template <typename TType>
pcrsdata<TType>::pcrsdata(std::vector<std::vector<TType> > &sizes, const std::vector<std::vector<TType> > &data)
: _crs(push_to_last(sizes)), _data(data)
{
	_iterator = new crsiterator<TType>(0, _crs.size() - 1, _crs.data(), _data.data());
	_titerators = new crsiterator<TType>*[threads()];
	#pragma omp parallel for
	for (size_t t = 0; t < _data.threads() - 1; t++) {
		_titerators[t] = new crsiterator<TType>(_crs.distribution()[t], _crs.distribution()[t + 1], _crs.data(), _data.data());
	}
	size_t t = _data.threads() - 1;
	_titerators[t] = new crsiterator<TType>(_crs.distribution()[t], _crs.distribution()[t + 1] - 1, _crs.data(), _data.data());

	std::vector<size_t> offsets(_crs.threads());

	#pragma omp parallel for
	for (size_t t = 0; t < _crs.threads(); t++) {
		size_t offset = 0;
		for (size_t i = _crs.distribution()[t]; i < _crs.distribution()[t + 1]; i++) {
			offset += _crs[i];
		}
		offsets[t] = offset;
	}

	size_t sum = 0;
	for (size_t i = 0; i < offsets.size(); i++) {
		size_t tmp = offsets[i];
		offsets[i] = sum;
		sum += tmp;
	}

	#pragma omp parallel for
	for (size_t t = 0; t < _crs.threads(); t++) {
		size_t offset = offsets[t];
		for (size_t i = _crs.distribution()[t]; i < _crs.distribution()[t + 1]; i++) {
			offset += _crs[i];
			_crs[i] = offset - _crs[i];
		}
	}

	sizes.back().pop_back();
}

template <typename TType>
pcrsdata<TType>::~pcrsdata()
{
	delete _iterator;
	for (size_t t = 0; t < threads(); t++) {
		delete _titerators[t];
	}
	delete[] _titerators;
}

}



#endif /* SRC_BASIS_CONTAINERS_PCRSDATA_HPP_ */
