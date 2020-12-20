
#include "matrix.csr.h"
#include "matrix.dense.h"
#include "vector.dense.h"
#include "vector.dense.distributed.h"
#include "math.h"
#include "esinfo/eslog.h"

#include <cstddef>

using namespace espreso;

MatrixCSR::MatrixCSR()
: _math(NULL)
{
	structureUpdated();
}

MatrixCSR::MatrixCSR(const MatrixCSR &other)
: Matrix(other), _math(NULL)
{
	deepCopy(&other);
}

MatrixCSR::MatrixCSR(const MatrixDense &other)
: Matrix(other), _math(NULL)
{
	*this = other;
}

MatrixCSR::MatrixCSR(esint nrows, esint ncols, esint nnz)
: DataMatrixCSR(nrows, ncols, nnz), _math(NULL)
{
	structureUpdated();
}

MatrixCSR& MatrixCSR::operator=(const MatrixCSR &other)
{
	if (this != &other) {
		deepCopy(&other);
	}
	return *this;
}

MatrixCSR& MatrixCSR::operator=(const MatrixDense &other)
{
	Matrix::_assign(&other);

	double eps = 1e-12;
	resize(other.nrows, other.ncols, other.nnz(eps));

	rows[0] = DataMatrixCSR::indexing;
	for (esint r = 0, i = 0; r < nrows; r++) {
		for (esint c = 0; c < ncols; c++) {
			if (other[r][c] < -eps || eps < other[r][c]) {
				cols[i] = c + DataMatrixCSR::indexing;
				vals[i] = other[r][c];
				++i;
			}
		}
		rows[r + 1] = i + DataMatrixCSR::indexing;
	}
	structureUpdated();
	return *this;
}

MatrixCSR::~MatrixCSR()
{
	if (_math) {
		delete _math;
	}
}

MatrixCSR* MatrixCSR::copy()
{
	return new MatrixCSR();
}

void MatrixCSR::structureUpdated()
{
	if (_math) { delete _math; }
	_math = new MATH::CSRHandler(nrows, ncols, nnz, rows, cols, vals);
}

void MatrixCSR::swap(Matrix *other)
{
	MatrixCSR *_other = other->downcast<MatrixCSR>();
	Matrix::_swap(other);
	DataMatrixCSR::swap(_other);
	MATH::CSRHandler *_handler = _math;
	_math = _other->_math;
	_other->_math = _handler;
}

void MatrixCSR::shallowCopy(const Matrix *other)
{
	Matrix::_assign(other);
	DataMatrixCSR::shallowCopy(other->downcast<MatrixCSR>());
	structureUpdated();
}

void MatrixCSR::shallowCopyStructure(const Matrix *other)
{
	Matrix::_assign(other);
	DataMatrixCSR::shallowCopyStructure(other->downcast<MatrixCSR>());
	structureUpdated();
}

void MatrixCSR::deepCopy(const Matrix *other)
{
	Matrix::_assign(other);
	DataMatrixCSR::deepCopy(other->downcast<MatrixCSR>());
	structureUpdated();
}

void MatrixCSR::deepCopyStructure(const Matrix *other)
{
	Matrix::_assign(other);
	DataMatrixCSR::deepCopyStructure(other->downcast<MatrixCSR>());
	structureUpdated();
}

void MatrixCSR::uniformCombination(const Matrix *first, const Matrix *second, int nfirst, int nsecond)
{
	DataMatrixCSR::uniformCombination(first->downcast<MatrixCSR>(), second->downcast<MatrixCSR>(), nfirst, nsecond);
	structureUpdated();
}

void MatrixCSR::fill(double value)
{
	DataMatrixCSR::fill(value);
}

void MatrixCSR::fillData(const Matrix *in)
{
	const MatrixCSR *_in = in->downcast<MatrixCSR>();
	DataMatrixCSR::fillValues(_in->nnz, _in->vals);
}

void MatrixCSR::fillCombinedData(const Matrix *in, esint roffset, esint coffset, esint nsize, esint sumsize)
{
	DataMatrixCSR::fillCombinedValues(in->downcast<MatrixCSR>(), roffset, coffset, nsize, sumsize);
}

void MatrixCSR::apply(const Vector *in, Vector *out)
{
	double *_in = NULL, *_out = NULL;
	if (dynamic_cast<const VectorDense*>(in)) { _in = dynamic_cast<const VectorDense*>(in)->vals; }
	if (dynamic_cast<const VectorDense*>(out)) { _out = dynamic_cast<const VectorDense*>(out)->vals; }
	if (dynamic_cast<const VectorDenseDistributed*>(in)) { _in = dynamic_cast<const VectorDenseDistributed*>(in)->vals; }
	if (dynamic_cast<const VectorDenseDistributed*>(out)) { _out = dynamic_cast<const VectorDenseDistributed*>(out)->vals; }

	if (_in == NULL || _out == NULL) {
		eslog::internalFailure("unsupported math operation.\n");
	}

	switch (type) {
	case MatrixType::REAL_SYMMETRIC_INDEFINITE:
	case MatrixType::REAL_SYMMETRIC_POSITIVE_DEFINITE:
		MATH::upCSRMatVecProduct(_math, _in, _out);
		break;
	case MatrixType::REAL_UNSYMMETRIC:
		MATH::CSRMatVecProduct(_math, _in, _out);
		break;
	}
}

void MatrixCSR::apply(const Vectors *in, Vectors *out)
{
	for (esint n = 0; n < in->nvectors; ++n) {
		apply(in->at(n), out->at(n));
	}
}

void MatrixCSR::scale(double alpha)
{
	MATH::vecScale(nnz, alpha, vals);
}

void MatrixCSR::add(double alpha, const Matrix *a)
{
	const MatrixCSR *_a = a->downcast<MatrixCSR>();
	MATH::vecAdd(nnz, vals, alpha, _a->vals);
}

void MatrixCSR::sum(double alpha, const Matrix *a, double beta, const Matrix *b)
{
	const MatrixCSR *_a = a->downcast<MatrixCSR>();
	const MatrixCSR *_b = b->downcast<MatrixCSR>();
	MATH::vecSum(nnz, vals, alpha, _a->vals, beta, _b->vals);
}

void MatrixCSR::addToCombination(double scale, const Matrix *in, esint roffset, esint coffset, esint nsize, esint sumsize)
{
	DataMatrixCSR::addToCombination(scale, in->downcast<MatrixCSR>(), roffset, coffset, nsize, sumsize);
}

void MatrixCSR::transpose()
{
	MatrixCSR tran(ncols, nrows, nnz);
	MATH::CSRTranspose(nrows, ncols, rows, cols, vals, tran.rows, tran.cols, tran.vals);
	tran.structureUpdated();
	tran.swap(this);
}

void MatrixCSR::transposeTo(MatrixCSR *out)
{
	out->resize(ncols, nrows, nnz);
	MATH::CSRTranspose(nrows, ncols, rows, cols, vals, out->rows, out->cols, out->vals);
	out->structureUpdated();
}

void MatrixCSR::fillDiagonal(Vector *diagonal) const
{
	DataMatrixCSR::fillDiagonal(diagonal->downcast<VectorDense>());
}

double MatrixCSR::norm()
{
	return MATH::vecNorm(nnz, vals);
}

void MatrixCSR::multiply(MatrixCSR &A, MatrixCSR &B, bool transposeA)
{
	MATH::CSRMatCSRMatProduct(_math, A._math, B._math, transposeA);
	DataMatrixCSR::operator=(_math);
}

void MatrixCSR::factorizeSymbolic()
{
	MATH::CSRMatFactorizeSymbolic(type, _math);
}

void MatrixCSR::factorizeNumeric()
{
	MATH::CSRMatFactorizeNumeric(type, _math);
}

void MatrixCSR::factorize()
{
	MATH::CSRMatFactorize(type, _math);
}

void MatrixCSR::solve(const MatrixDense &rhs, MatrixDense &solution)
{
	MATH::CSRMatSolve(type, _math, rhs.nrows, rhs.vals, solution.vals);
}

void MatrixCSR::solve(const VectorDense &rhs, VectorDense &solution)
{
	MATH::CSRMatSolve(type, _math, 1, rhs.vals, solution.vals);
}

void MatrixCSR::removeLower(MatrixType type)
{
	this->type = type;
	allowUpdating();
	MATH::CSRRemoveLower(nrows, ncols, rows, cols, vals);
	nnz = rows[nrows] - rows[0];
	structureUpdated();
}
