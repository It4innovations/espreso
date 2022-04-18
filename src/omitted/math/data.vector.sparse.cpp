
#include "data.vector.sparse.h"
#include "data.vector.dense.h"
#include "math.h"

#include <cstddef>
#include <cstring>

using namespace espreso;

_DataVectorSparse::_DataVectorSparse()
: size(0), nnz(0), indices(NULL), vals(NULL)
{

}
void _DataVectorSparse::alloc(esint size, esint nnz)
{
	clear();
	this->size = size;
	this->nnz = nnz;

	if (nnz) {
		indices = new esint[nnz];
		vals = new double[nnz];
	}
}

void _DataVectorSparse::clear()
{
	size = 0;
	nnz = 0;
	if (indices) { delete[] indices; indices = NULL; }
	if (vals) { delete[] vals; vals = NULL; }
}

void DataVectorSparse::combineIndices(esint *result, esint *first, esint *second, esint *firstend, esint *secondend, esint nfirst, esint nsecond)
{
	int sum = nfirst + nsecond;
	esint *i1 = first;
	esint *i2 = second;
	esint *i = result;
	while (i1 != firstend || i2 != secondend) {
		if (i1 == firstend) {
			*i = sum * (*i2 / nsecond) + *i2 % nsecond + nfirst;
			++i; ++i2;
			continue;
		}
		if (i2 == secondend) {
			*i = sum * (*i1 / nfirst) + *i1 % nfirst;
			++i; ++i1;
			continue;
		}
		if (*i1 / nfirst <= *i2 / nsecond) {
			*i = sum * (*i1 / nfirst) + *i1 % nfirst;
			++i; ++i1;
		} else {
			*i = sum * (*i2 / nsecond) + *i2 % nsecond + nfirst;
			++i; ++i2;
		}
	}
}

DataVectorSparse::DataVectorSparse()
{

}

DataVectorSparse::DataVectorSparse(esint size, esint nnz)
{
	resize(size, nnz);
}

DataVectorSparse::~DataVectorSparse()
{
	_allocated.clear();
}

void DataVectorSparse::shiftData(esint offset)
{
	vals += offset;
}

void DataVectorSparse::resize(esint size, esint nnz)
{
	_allocated.alloc(size, nnz);
	_DataVectorSparse::operator=(_allocated);
}

void DataVectorSparse::swap(DataVectorSparse *other)
{
	_DataVectorSparse _data = _allocated;
	_allocated = other->_allocated;
	other->_allocated = _data;

	_data = *this;
	this->_DataVectorSparse::operator=(*other);
	other->_DataVectorSparse::operator=(_data);
}

void DataVectorSparse::shallowCopy(const DataVectorSparse *other)
{
	_allocated.clear();
	_DataVectorSparse::operator=(*other);
}

void DataVectorSparse::shallowCopyStructure(const DataVectorSparse *other)
{
	shallowCopy(other);
	_allocated.vals = vals = new double[nnz];
}

void DataVectorSparse::shallowCopyFromHolder(const DataVectorSparse *other, esint offset, esint nvectors)
{
	shallowCopy(other);
	size /= nvectors;
	nnz /= nvectors;
	shiftData(offset * nnz);
}

void DataVectorSparse::deepCopy(const DataVectorSparse *other)
{
	deepCopyStructure(other);
	fillValues(other->vals);
}

void DataVectorSparse::deepCopyStructure(const DataVectorSparse *other)
{
	_allocated.clear();
	_allocated.alloc(other->size, other->nnz);
	_DataVectorSparse::operator=(_allocated);
	fillPattern(other->indices);
}

void DataVectorSparse::uniformCombination(const DataVectorSparse *first, const DataVectorSparse *second, int nfirst, int nsecond)
{
	resize(first->size + second->size, first->nnz + second->nnz);
	combineIndices(indices, first->indices, second->indices, first->indices + first->nnz, second->indices + second->nnz, nfirst, nsecond);
}

void DataVectorSparse::fill(double value)
{
	for (esint i = 0; i < nnz; ++i) {
		vals[i] = value;
	}
}

void DataVectorSparse::fillPattern(esint *indices)
{
	memcpy(this->indices, indices, sizeof(esint) * nnz);
}

void DataVectorSparse::fillValues(double *vals)
{
	memcpy(this->vals, vals, sizeof(double) * nnz);
}

void DataVectorSparse::fillDenseValues(double *vals)
{
	for (esint n = 0 ; n < nnz; ++n) {
		this->vals[n] = vals[indices[n]];
	}
}

void DataVectorSparse::fillCombinedValues(const DataVectorSparse *in, esint offset, esint nsize, esint sumsize)
{
	// TODO: generalize
	for (esint i = 0; i < in->nnz / nsize; ++i) {
		for (esint j = 0; j < nsize; ++j) {
			vals[i * sumsize + offset + j] = in->vals[i * nsize + j];
		}
	}
}

void DataVectorSparse::fillValuesFromCombination(const DataVectorSparse *in, esint offset, esint nsize, esint sumsize)
{
	// TODO: generalize
	for (esint i = 0; i < nnz / nsize; ++i) {
		for (esint j = 0; j < nsize; ++j) {
			vals[i * nsize + j] = in->vals[i * sumsize + offset + j];
		}
	}
}

void DataVectorSparse::addToCombination(double alpha, const DataVectorSparse *in, esint offset, esint nsize, esint sumsize)
{
	// TODO
}

