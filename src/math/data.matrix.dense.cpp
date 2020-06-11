
#include "data.matrix.dense.h"
#include "data.vector.dense.h"
#include "math.h"

#include <cstddef>
#include <cstring>

using namespace espreso;

_DataMatrixDense::_DataMatrixDense()
: nrows(0), ncols(0), vals(NULL), _maxvals(0)
{

}

void _DataMatrixDense::realloc(esint nrows, esint ncols)
{
	if (_maxvals < nrows * ncols) {
		clear();
		_maxvals = nrows * ncols;
		vals = new double[_maxvals];
	}
	this->nrows = nrows;
	this->ncols = ncols;
}

void _DataMatrixDense::clear()
{
	nrows = ncols = _maxvals = 0;
	if (vals) { delete[] vals; vals = NULL; }
}

DataMatrixDense::DataMatrixDense()
{

}

DataMatrixDense::DataMatrixDense(esint nrows, esint ncols)
{
	resize(nrows, ncols);
	fill(0); // TODO: remove (change kernels)
}

DataMatrixDense::DataMatrixDense(esint nrows, esint ncols, double *vals)
{
	set(nrows, ncols, vals);
}

DataMatrixDense::DataMatrixDense(const DataMatrixDense &other)
{
	deepCopy(&other);
}

DataMatrixDense::~DataMatrixDense()
{
	_allocated.clear();
}

void DataMatrixDense::set(esint nrows, esint ncols, double *vals)
{
	_allocated.clear();
	this->nrows = nrows;
	this->ncols = ncols;
	this->vals = vals;
}

void DataMatrixDense::resize(esint nrows, esint ncols)
{
	_allocated.realloc(nrows, ncols);
	_DataMatrixDense::operator=(_allocated);
}

void DataMatrixDense::swap(DataMatrixDense *other)
{
	_DataMatrixDense _data = _allocated;
	_allocated = other->_allocated;
	other->_allocated = _data;

	_data = *this;
	this->_DataMatrixDense::operator=(*other);
	other->_DataMatrixDense::operator=(_data);
}

void DataMatrixDense::shallowCopy(const DataMatrixDense *other)
{
	_allocated.clear();
	_DataMatrixDense::operator=(*other);
}

void DataMatrixDense::deepCopy(const DataMatrixDense *other)
{
	_allocated.realloc(other->nrows, other->ncols);
	_DataMatrixDense::operator=(_allocated);
	fillValues(other->vals);
}

void DataMatrixDense::uniformCombination(const DataMatrixDense *first, const DataMatrixDense *second, int nfirst, int nsecond)
{
	nrows = first->nrows + second->nrows;
	ncols = first->ncols + second->ncols;
	resize(nrows, ncols);
}

void DataMatrixDense::fill(double value)
{
	for (esint i = 0; i < nrows * ncols; ++i) {
		vals[i] = value;
	}
}

void DataMatrixDense::fillValues(double *vals)
{
	memcpy(this->vals, vals, sizeof(double) * nrows * ncols);
}

void DataMatrixDense::fillCombinedValues(const DataMatrixDense *in,  esint roffset, esint coffset, esint rsize, esint csize, esint rsum, esint csum)
{
	for (esint r = 0; r < in->nrows / rsize; ++r) {
		for (esint i = 0; i < rsize; ++i) {
			for (esint c = 0; c < in->ncols / csize; ++c) {
				for (esint j = 0; j < csize; ++j) {
					vals[ncols * (r * rsum + roffset + i) + c * csum + coffset + j] = in->vals[in->ncols * (r * rsize + i) + c * csize + j];
				}
			}
		}
	}
}

void DataMatrixDense::addToCombination(double alpha, const DataMatrixDense *in,  esint roffset, esint coffset, esint rsize, esint csize, esint rsum, esint csum)
{
	for (esint r = 0; r < in->nrows / rsize; ++r) {
		for (esint i = 0; i < rsize; ++i) {
			for (esint c = 0; c < in->ncols / csize; ++c) {
				for (esint j = 0; j < csize; ++j) {
					vals[ncols * (r * rsum + roffset + i) + c * rsum + coffset + j] += in->vals[in->ncols * (r * rsize + i) + c * csize + j];
				}
			}
		}
	}
}

void DataMatrixDense::fillDiagonal(DataVectorDense *diagonal) const
{
	for (esint r = 0; r < nrows; ++r) {
		diagonal[r] = vals[r * ncols + r];
	}
}



