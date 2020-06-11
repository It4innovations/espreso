
#include "matrix.ijv.h"
#include "matrix.dense.h"
#include "vector.dense.h"
#include "math.h"
#include "esinfo/eslog.h"

#include <cstddef>

using namespace espreso;

MatrixIJV::MatrixIJV()
: _math(NULL)
{
	structureUpdated();
}

MatrixIJV::MatrixIJV(const MatrixIJV &other)
: Matrix(other), _math(NULL)
{
	deepCopy(&other);
}

MatrixIJV::MatrixIJV(const MatrixDense &other)
: Matrix(other), _math(NULL)
{
	*this = other;
}

MatrixIJV::MatrixIJV(esint nrows, esint ncols, esint nnz)
: DataMatrixIJV(nrows, ncols, nnz), _math(NULL)
{
	structureUpdated();
}

MatrixIJV& MatrixIJV::operator=(const MatrixIJV &other)
{
	if (this != &other) {
		deepCopy(&other);
	}
	return *this;
}

MatrixIJV& MatrixIJV::operator=(const MatrixDense &other)
{
	Matrix::_assign(&other);

	double eps = 1e-12;
	resize(other.nrows, other.ncols, other.nnz(eps));

	for (esint r = 0, i = 0; r < nrows; r++) {
		for (esint c = 0; c < ncols; c++) {
			if (other[r][c] < -eps || eps < other[r][c]) {
				cols[i] = r + DataMatrixIJV::indexing;
				cols[i] = c + DataMatrixIJV::indexing;
				vals[i] = other[r][c];
				++i;
			}
		}
	}

	return *this;
}

MatrixIJV::~MatrixIJV()
{
	if (_math) { delete _math; }
}

MatrixIJV* MatrixIJV::copy()
{
	return new MatrixIJV();
}

void MatrixIJV::structureUpdated()
{
	if (_math) { delete _math; }
	_math = new MATH::IJVHandler(nrows, ncols, nnz, rows, cols, vals);
}

void MatrixIJV::swap(Matrix *other)
{
	MatrixIJV *_other = other->downcast<MatrixIJV>();
	Matrix::_swap(other);
	DataMatrixIJV::swap(_other);
	MATH::IJVHandler *_handler = _math;
	_math = _other->_math;
	_other->_math = _handler;
}

void MatrixIJV::shallowCopy(const Matrix *other)
{
	Matrix::_assign(other);
	DataMatrixIJV::shallowCopy(other->downcast<MatrixIJV>());
	structureUpdated();
}

void MatrixIJV::shallowCopyStructure(const Matrix *other)
{
	Matrix::_assign(other);
	DataMatrixIJV::shallowCopyStructure(other->downcast<MatrixIJV>());
	structureUpdated();
}

void MatrixIJV::deepCopy(const Matrix *other)
{
	Matrix::_assign(other);
	DataMatrixIJV::deepCopy(other->downcast<MatrixIJV>());
	structureUpdated();
}

void MatrixIJV::deepCopyStructure(const Matrix *other)
{
	Matrix::_assign(other);
	DataMatrixIJV::deepCopyStructure(other->downcast<MatrixIJV>());
	structureUpdated();
}

void MatrixIJV::uniformCombination(const Matrix *first, const Matrix *second, int nfirst, int nsecond)
{
	DataMatrixIJV::uniformCombination(first->downcast<MatrixIJV>(), second->downcast<MatrixIJV>(), nfirst, nsecond);
	structureUpdated();
}

void MatrixIJV::fill(double value)
{
	DataMatrixIJV::fill(value);
}

void MatrixIJV::fillData(const Matrix *in)
{
	const MatrixIJV *_in = in->downcast<MatrixIJV>();
	DataMatrixIJV::fillValues(_in->nnz, _in->vals);
}

void MatrixIJV::fillCombinedData(const Matrix *in,  esint roffset, esint coffset, esint rsize, esint csize, esint rsum, esint csum)
{
	DataMatrixIJV::fillCombinedValues(in->downcast<MatrixIJV>(), roffset, coffset, rsize, csize, rsum, csum);
}

void MatrixIJV::apply(const Vector *in, Vector *out)
{

}

void MatrixIJV::apply(const Vectors *in, Vectors *out)
{

}

void MatrixIJV::scale(double alpha)
{
	MATH::vecScale(nnz, alpha, vals);
}

void MatrixIJV::add(double alpha, const Matrix *a)
{
	const MatrixIJV *_a = a->downcast<MatrixIJV>();
	MATH::vecAdd(nnz, vals, alpha, _a->vals);
}

void MatrixIJV::sum(double alpha, const Matrix *a, double beta, const Matrix *b)
{
	const MatrixIJV *_a = a->downcast<MatrixIJV>();
	const MatrixIJV *_b = b->downcast<MatrixIJV>();
	MATH::vecSum(nnz, vals, alpha, _a->vals, beta, _b->vals);
}

void MatrixIJV::addToCombination(double scale, const Matrix *in,  esint roffset, esint coffset, esint rsize, esint csize, esint rsum, esint csum)
{
	DataMatrixIJV::addToCombination(scale, in->downcast<MatrixIJV>(), roffset, coffset, rsize, csize, rsum, csum);
}

void MatrixIJV::transpose()
{
	esint _t = nrows;
	nrows = ncols;
	ncols = _t;
	esint *d = rows;
	rows = cols;
	cols = d;
}

void MatrixIJV::transposeTo(MatrixIJV *out)
{
	out->DataMatrixIJV::resize(ncols, nrows, nnz);
	out->DataMatrixIJV::fillPattern(nnz, cols, rows);
	out->DataMatrixIJV::fillValues(nnz, vals);
}

void MatrixIJV::fillDiagonal(Vector *diagonal) const
{
	DataMatrixIJV::fillDiagonal(diagonal->downcast<VectorDense>());
}

double MatrixIJV::norm()
{
	return MATH::vecNorm(nnz, vals);
}

void MatrixIJV::removeLower(MatrixType type)
{
	eslog::internalFailure("call empty function.\n");
}

