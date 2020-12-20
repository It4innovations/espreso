
#include "data.matrix.csr.h"
#include "data.vector.dense.h"
#include "math.h"
#include "esinfo/eslog.h"

#include <cstddef>
#include <cstring>
#include <vector>

using namespace espreso;

const esint DataMatrixCSR::indexing = 1;

_DataMatrixCSR::_DataMatrixCSR()
: nrows(0), ncols(0), nnz(0), rows(NULL), cols(NULL), vals(NULL)
{

}

void _DataMatrixCSR::alloc(esint nrows, esint ncols, esint nnz)
{
	clear();
	this->nrows = nrows;
	this->ncols = ncols;
	this->nnz = nnz;

	if (nnz) {
		rows = new esint[nrows + 1];
		cols = new esint[nnz];
		vals = new double[nnz];
	}
}

void _DataMatrixCSR::clear()
{
	nrows = 0;
	ncols = 0;
	nnz = 0;
	if (rows) { delete[] rows; rows = NULL; }
	if (cols) { delete[] cols; cols = NULL; }
	if (vals) { delete[] vals; vals = NULL; }
}

void DataMatrixCSR::combineIndices(esint *result, esint *first, esint *second, esint *firstend, esint *secondend, esint nfirst, esint nsecond)
{
	eslog::internalFailure("call empty function.\n");
}

DataMatrixCSR::DataMatrixCSR()
{

}

DataMatrixCSR::DataMatrixCSR(const DataMatrixCSR &other)
{
	deepCopy(&other);
}

DataMatrixCSR::DataMatrixCSR(const MATH::CSRHandler *handler)
{
	_allocated.clear();
	handler->sizes(nrows, ncols, nnz);
	handler->data(rows, cols, vals);
}

DataMatrixCSR::DataMatrixCSR(DataMatrixCSR &&other)
{
	shallowCopy(&other);
	_allocated = other._allocated;
	other._allocated = _DataMatrixCSR();
}

DataMatrixCSR::DataMatrixCSR(esint nrows, esint ncols, esint nnz)
{
	resize(nrows, ncols, nnz);
}

DataMatrixCSR& DataMatrixCSR::operator=(const DataMatrixCSR &other)
{
	if (this != &other) {
		deepCopy(&other);
	}
	return *this;
}

DataMatrixCSR& DataMatrixCSR::operator=(const MATH::CSRHandler *handler)
{
	_allocated.clear();
	handler->sizes(nrows, ncols, nnz);
	handler->data(rows, cols, vals);
	return *this;
}

DataMatrixCSR::~DataMatrixCSR()
{
	_allocated.clear();
}

void DataMatrixCSR::allowUpdating()
{
	if (_allocated.rows == NULL || _allocated.cols == NULL || _allocated.vals == NULL) {
		DataMatrixCSR other;
		other.deepCopy(this);
		swap(&other);
	}
}

void DataMatrixCSR::resize(esint nrows, esint ncols, esint nnz)
{
	_allocated.clear();
	_allocated.alloc(nrows, ncols, nnz);
	_DataMatrixCSR::operator=(_allocated);
}

void DataMatrixCSR::swap(DataMatrixCSR *other)
{
	_DataMatrixCSR _data = _allocated;
	_allocated = other->_allocated;
	other->_allocated = _data;

	_data = *this;
	this->_DataMatrixCSR::operator=(*other);
	other->_DataMatrixCSR::operator=(_data);
}

void DataMatrixCSR::shallowCopy(const DataMatrixCSR *other)
{
	_allocated.clear();
	_DataMatrixCSR::operator=(*other);
}

void DataMatrixCSR::shallowCopyStructure(const DataMatrixCSR *other)
{
	shallowCopy(other);
	_allocated.vals = vals = new double[nnz];
}

void DataMatrixCSR::deepCopy(const DataMatrixCSR *other)
{
	deepCopyStructure(other);
	fillValues(other->nnz, other->vals);
}

void DataMatrixCSR::deepCopyStructure(const DataMatrixCSR *other)
{
	_allocated.clear();
	_allocated.alloc(other->nrows, other->ncols, other->nnz);
	_DataMatrixCSR::operator=(_allocated);
	fillPattern(other->nrows, other->rows, other->cols);
}

void DataMatrixCSR::fillPattern(esint nrows, esint *rows, esint *cols)
{
	if (nrows) {
		memcpy(this->rows, rows, sizeof(esint) * (nrows + 1));
		memcpy(this->cols, cols, sizeof(esint) * (rows[nrows] - rows[0]));
	}
}

void DataMatrixCSR::fillValues(esint nnz, double *vals)
{
	if (nnz) {
		memcpy(this->vals, vals, sizeof(double) * nnz);
	}
}

void DataMatrixCSR::uniformCombination(const DataMatrixCSR *first, const DataMatrixCSR *second, int nfirst, int nsecond)
{
// TODO: generalize
	if (first->nrows % nfirst != 0 || second->nrows % nsecond != 0) {
		eslog::internalFailure("cannot combine matrices.");
	}

	int sum = nfirst + nsecond;

	std::vector<esint> newrows, newcols;
	newrows.push_back(1);
	for (esint n = 0; n < first->nrows / nfirst; ++n) {
		for (int r = 0; r < sum; ++r) {
			for (esint i = first->rows[n * nfirst]; i < first->rows[n * nfirst + 1]; i += nfirst) {
				esint column = first->cols[i - first->rows[0]] - first->rows[0];
				for (int c = 0; c < sum; c++) {
					newcols.push_back(sum * (column / nfirst) + c + newrows[0]);
				}
			}
			newrows.push_back(newcols.size() + newrows[0]);
		}
	}

	_allocated.clear();
	_allocated.alloc(first->nrows + second->nrows, first->ncols + second->ncols, newcols.size());
	_DataMatrixCSR::operator=(_allocated);

	if (first->nnz + second->nnz) {
		memcpy(rows, newrows.data(), newrows.size() * sizeof(esint));
		memcpy(cols, newcols.data(), newcols.size() * sizeof(esint));
	}
}

void DataMatrixCSR::fill(double value)
{
	for (esint i = 0; i < nnz; ++i) {
		vals[i] = value;
	}
}

void DataMatrixCSR::fillCombinedValues(const DataMatrixCSR *in, esint roffset, esint coffset, esint nsize, esint sumsize)
{
	eslog::internalFailure("call empty function.\n");
}

void DataMatrixCSR::addToCombination(double scale, const DataMatrixCSR *in, esint roffset, esint coffset, esint nsize, esint sumsize)
{
	for (esint r = 0; r < in->nrows; r++) {
		for (esint c = 0; c < in->rows[r + 1] - in->rows[r]; c++) {
			esint tr = sumsize * (r / nsize) + roffset + r % nsize;
			esint rc = sumsize * (c / nsize) + coffset + c % nsize;
			vals[rows[tr] + rc - rows[0]] += scale * in->vals[in->rows[r] + c - in->rows[0]];
		}
	}
}

void DataMatrixCSR::fillDiagonal(DataVectorDense *diagonal) const
{
	for (esint r = 0; r < nrows; r++) {
		for (esint c = rows[r] - 1; c < rows[r + 1] - 1; c++) {
			if (cols[c] == r + 1) {
				diagonal->vals[r] = vals[c];
				break;
			}
		}
	}
}

