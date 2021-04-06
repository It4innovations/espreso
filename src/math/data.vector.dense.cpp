
#include "data.vector.dense.h"
#include "data.vector.sparse.h"
#include "math.h"

#include <cstddef>
#include <cstring>

using namespace espreso;

_DataVectorDense::_DataVectorDense()
: size(0), vals(NULL), _maxsize(0)
{

}

void _DataVectorDense::realloc(esint size)
{
	if (_maxsize < size) {
		clear();
		_maxsize = size;
		vals = new double[_maxsize];
	}
	this->size = size;
}

void _DataVectorDense::clear()
{
	size = _maxsize = 0;
	if (vals) { delete[] vals; vals = NULL; }
}

DataVectorDense::DataVectorDense()
{

}

DataVectorDense::DataVectorDense(esint nvals)
{
	resize(nvals);
}

DataVectorDense::DataVectorDense(esint size, double *vals)
{
	this->size = size;
	this->vals = vals;
}

DataVectorDense::~DataVectorDense()
{
	_allocated.clear();
}

void DataVectorDense::shiftData(esint offset)
{
	vals += offset;
}

void DataVectorDense::resize(esint size)
{
	_allocated.realloc(size);
	_DataVectorDense::operator=(_allocated);
}

void DataVectorDense::swap(DataVectorDense *other)
{
	_DataVectorDense _data = _allocated;
	_allocated = other->_allocated;
	other->_allocated = _data;

	_data = *this;
	this->_DataVectorDense::operator=(*other);
	other->_DataVectorDense::operator=(_data);
}

void DataVectorDense::shallowCopy(const DataVectorDense *other)
{
	_allocated.clear();
	_DataVectorDense::operator=(*other);
}

void DataVectorDense::shallowCopyStructure(const DataVectorDense *other)
{
	_allocated.clear();
	_DataVectorDense::operator=(*other);
	_allocated.vals = vals = new double[other->size];
}

void DataVectorDense::shallowCopyFromHolder(const DataVectorDense *other, esint offset, esint nvectors)
{
	shallowCopy(other);
	size /= nvectors;
	shiftData(offset * size);
}

void DataVectorDense::deepCopy(const DataVectorDense *other)
{
	deepCopyStructure(other);
	fillValues(other->vals);
}

void DataVectorDense::deepCopyStructure(const DataVectorDense *other)
{
	_allocated.realloc(other->size);
	_DataVectorDense::operator=(_allocated);
}

void DataVectorDense::uniformCombination(const DataVectorDense *first, const DataVectorDense *second)
{
	size = first->size + second->size;
	resize(size);
}

void DataVectorDense::fill(double value)
{
	for (esint i = 0; i < size; ++i) {
		vals[i] = value;
	}
}

void DataVectorDense::fillValues(double *vals)
{
	memcpy(this->vals, vals, sizeof(double) * size);
}

void DataVectorDense::fillSparseValues(esint nnz, esint *indices, double *vals)
{
	for (esint i = 0; i < nnz; i++) {
		this->vals[indices[i]] = vals[i];
	}
}

void DataVectorDense::fillCombinedValues(const DataVectorDense *in, esint offset, esint nsize, esint sumsize)
{
	for (esint i = 0; i < in->size / nsize; ++i) {
		for (esint j = 0; j < nsize; ++j) {
			vals[i * sumsize + offset + j] = in->vals[i * nsize + j];
		}
	}
}

void DataVectorDense::fillValuesFromCombination(const DataVectorDense *in, esint offset, esint nsize, esint sumsize)
{
	for (esint i = 0; i < size / nsize; ++i) {
		for (esint j = 0; j < nsize; ++j) {
			vals[i * nsize + j] = in->vals[i * sumsize + offset + j];
		}
	}
}

void DataVectorDense::addToCombination(double alpha, const DataVectorDense *in, esint offset, esint nsize, esint sumsize)
{
	for (esint i = 0; i < in->size / nsize; ++i) {
		for (esint j = 0; j < nsize; ++j) {
			vals[i * sumsize + offset + j] += alpha * in->vals[i * nsize + j];
		}
	}
}

