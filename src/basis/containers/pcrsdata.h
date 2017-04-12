
#ifndef SRC_BASIS_CONTAINERS_PCRSDATA_H_
#define SRC_BASIS_CONTAINERS_PCRSDATA_H_

#include <vector>

#include "parray.h"
#include "crsiterator.h"

/**
 * Parallel Compressed Row Storage format class.
 */

namespace espreso {

template <typename TType>
class pcrsdata {

public:
	pcrsdata(parray<TType> &&crs, parray<TType> &&data);
	pcrsdata(std::vector<std::vector<TType> > &sizes, const std::vector<std::vector<TType> > &data);

	size_t size() const { return _data.size(); }
	size_t size(size_t index) { return _crs[index + 1] - _crs[index]; }
	size_t threads() const { return _data.threads(); }

	crsiterator<TType>& iterator() { return *_iterator; }
	crsiterator<TType>& iterator(size_t thread) { return *_titerators[thread]; }
	TType* crs() { return _crs.data(); }
	TType* data() { return _data.data(); }
	TType* data(size_t index) { return _data.data() + _crs[index]; }
	parray<TType>& pcrs() { return _crs; }
	parray<TType>& pdata() { return _data; }

	const crsiterator<TType>& iterator() const { return *_iterator; }
	const crsiterator<TType>& iterator(size_t thread) const { return *_titerators[thread]; }
	const TType* crs() const { return _crs.data(); }
	const TType* data() const { return _data.data(); }
	const TType* data(size_t index) const { return _data.data() + _crs[index]; }
	const parray<TType>& pcrs() const { return _crs; }
	const parray<TType>& pdata() const { return _data; }

	~pcrsdata();

private:
	crsiterator<TType>* _iterator;
	crsiterator<TType>** _titerators;

	parray<TType> _crs;
	parray<TType> _data;

};

}

#include "pcrsdata.hpp"


#endif /* SRC_BASIS_CONTAINERS_PCRSDATA_H_ */
