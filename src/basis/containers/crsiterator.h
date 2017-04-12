
#ifndef SRC_BASIS_CONTAINERS_CRSITERATOR_H_
#define SRC_BASIS_CONTAINERS_CRSITERATOR_H_

#include <cstddef>

/**
 * Wrapper for data stored in CSR format.
 *
 * It mimics std::vector methods that do not modify the data.
 * The instance of this class always points to a row.
 * By the next method, the user is able to move to the next row.
 * When next returns false, the class points again to the first row.
 */

namespace espreso {

template <typename TType>
class crsiterator {

public:
	crsiterator(size_t begin, size_t end, TType *offsets, const TType *data)
	: _origin(offsets), _begin(offsets + begin), _current(offsets + end), _end(offsets + end), _data(data) {}

	bool next() const { return _end != (_current = _current == _end ? _begin : _current + 1); }
	void seek(size_t iterator) const { _current = _origin + iterator; }
	void reset() const { _current = _end; }

	size_t size() const { return *(_current + 1) - *_current; }
	size_t index() const { return _current - _origin; }

	TType& operator[] (size_t n) { return _data[*_current + n]; }
	TType* data() { return _data + *_current; }
	TType& back() { return _data[*(_current + 1) - 1]; }
	TType& front() { return _data[*_current]; }
	TType* begin() { return _data + *_current; }
	TType* end() { return _data + *(_current + 1); }

	const TType& operator[] (size_t n) const { return _data[*_current + n]; }
	const TType* data() const { return _data + *_current; }
	const TType& back() const { return _data[*(_current + 1) - 1]; }
	const TType& front() const { return _data[*_current]; }
	const TType* begin() const { return _data + *_current; }
	const TType* end() const { return _data + *(_current + 1); }
	const TType* cbegin() const { return _data + *_current; }
	const TType* cend() const { return _data + *(_current + 1); }

private:
	TType *_origin;
	TType *_begin;
	mutable TType *_current;
	TType *_end;

	const TType *_data;
};

}



#endif /* SRC_BASIS_CONTAINERS_CRSITERATOR_H_ */
