
#include "data.matrix.ijv.h"
#include "data.vector.dense.h"
#include "math.h"
#include "esinfo/eslog.h"

#include <cstddef>
#include <cstring>
#include <vector>

using namespace espreso;

const esint DataMatrixIJV::indexing = 1;

_DataMatrixIJV::_DataMatrixIJV()
: nrows(0), ncols(0), nnz(0), rows(NULL), cols(NULL), vals(NULL)
{

}

void _DataMatrixIJV::alloc(esint nrows, esint ncols, esint nnz)
{
	clear();
	this->nrows = nrows;
	this->ncols = ncols;
	this->nnz = nnz;

	if (nnz) {
		rows = new esint[nnz];
		cols = new esint[nnz];
		vals = new double[nnz];
	}
}

void _DataMatrixIJV::clear()
{
	nrows = 0;
	ncols = 0;
	nnz = 0;
	if (rows) { delete[] rows; rows = NULL; }
	if (cols) { delete[] cols; cols = NULL; }
	if (vals) { delete[] vals; vals = NULL; }
}

void DataMatrixIJV::combineIndices(esint *result, esint *first, esint *second, esint *firstend, esint *secondend, esint nfirst, esint nsecond)
{
	eslog::internalFailure("call empty function.\n");
}

DataMatrixIJV::DataMatrixIJV()
{

}

DataMatrixIJV::DataMatrixIJV(const DataMatrixIJV &other)
{
	deepCopy(&other);
}

DataMatrixIJV::DataMatrixIJV(DataMatrixIJV &&other)
{
	shallowCopy(&other);
	_allocated = other._allocated;
	other._allocated = _DataMatrixIJV();
}

DataMatrixIJV::DataMatrixIJV(esint nrows, esint ncols, esint nnz)
{
	resize(nrows, ncols, nnz);
}

DataMatrixIJV& DataMatrixIJV::operator=(const DataMatrixIJV &other)
{
	if (this != &other) {
		deepCopy(&other);
	}
	return *this;
}

DataMatrixIJV::~DataMatrixIJV()
{
	_allocated.clear();
}

void DataMatrixIJV::swap(DataMatrixIJV *other)
{
	_DataMatrixIJV _data = _allocated;
	_allocated = other->_allocated;
	other->_allocated = _data;

	_data = *this;
	this->_DataMatrixIJV::operator=(*other);
	other->_DataMatrixIJV::operator=(_data);
}

void DataMatrixIJV::shallowCopy(const DataMatrixIJV *other)
{
	_allocated.clear();
	_DataMatrixIJV::operator=(*other);
}

void DataMatrixIJV::shallowCopyStructure(const DataMatrixIJV *other)
{
	shallowCopy(other);
	_allocated.vals = vals = new double[nnz];
}

void DataMatrixIJV::deepCopy(const DataMatrixIJV *other)
{
	deepCopyStructure(other);
	fillValues(other->nnz, other->vals);
}

void DataMatrixIJV::deepCopyStructure(const DataMatrixIJV *other)
{
	_allocated.clear();
	_allocated.alloc(other->nrows, other->ncols, other->nnz);
	_DataMatrixIJV::operator=(_allocated);
	fillPattern(other->nrows, other->rows, other->cols);
}

void DataMatrixIJV::fillPattern(esint nnz, esint *rows, esint *cols)
{
	if (nnz) {
		memcpy(this->rows, rows, sizeof(esint) * nnz);
		memcpy(this->cols, cols, sizeof(esint) * nnz);
	}
}

void DataMatrixIJV::fillValues(esint nnz, double *vals)
{
	if (nnz) {
		memcpy(this->vals, vals, sizeof(double) * nnz);
	}
}

void DataMatrixIJV::uniformCombination(const DataMatrixIJV *first, const DataMatrixIJV *second, int nfirst, int nsecond)
{
	eslog::internalFailure("call empty function.\n");
}

void DataMatrixIJV::allowUpdating()
{
	if (_allocated.rows == NULL || _allocated.cols == NULL || _allocated.vals == NULL) {
		DataMatrixIJV other;
		other.deepCopy(this);
		swap(&other);
	}
}

void DataMatrixIJV::resize(esint nrows, esint ncols, esint nnz)
{
	_allocated.clear();
	_allocated.alloc(nrows, ncols, nnz);
	_DataMatrixIJV::operator=(_allocated);
}

void DataMatrixIJV::fill(double value)
{
	for (esint i = 0; i < nnz; ++i) {
		vals[i] = value;
	}
}

void DataMatrixIJV::fillCombinedValues(const DataMatrixIJV *in,  esint roffset, esint coffset, esint rsize, esint csize, esint rsum, esint csum)
{
	eslog::internalFailure("call empty function.\n");
}

void DataMatrixIJV::addToCombination(double scale, const DataMatrixIJV *in,  esint roffset, esint coffset, esint rsize, esint csize, esint rsum, esint csum)
{
	eslog::internalFailure("call empty function.\n");
}

void DataMatrixIJV::fillDiagonal(DataVectorDense *diagonal) const
{
	eslog::internalFailure("call empty function.\n");
}


