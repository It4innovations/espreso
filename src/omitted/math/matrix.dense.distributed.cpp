
#include "matrix.dense.distributed.h"
#include "matrix.csr.distributed.h"
#include "matrix.dense.h"
#include "vector.dense.h"
#include "vector.dense.distributed.h"
#include "data.synchronization.h"
#include "data.mv.h"
#include "math.h"

#include "esinfo/eslog.hpp"
#include "esinfo/mpiinfo.h"

#include <cstddef>

using namespace espreso;

MatrixDenseDistributed::MatrixDenseDistributed()
{

}

MatrixDenseDistributed::MatrixDenseDistributed(esint nrows, esint ncols, esint nhalo, esint nneighbors)
: DataMatrixDense(nrows, ncols), DataDistributed(nhalo, info::mpi::size, nneighbors)
{
	structureUpdated();
}

MatrixDenseDistributed::MatrixDenseDistributed(const MatrixDenseDistributed &other)
{
	deepCopy(&other);
}

MatrixDenseDistributed& MatrixDenseDistributed::operator=(const MatrixDenseDistributed &other)
{
	if (this != &other) {
		deepCopy(&other);
	}
	return *this;
}

MatrixDenseDistributed::~MatrixDenseDistributed()
{

}


MatrixDenseDistributed* MatrixDenseDistributed::copy()
{
	return new MatrixDenseDistributed();;
}

void MatrixDenseDistributed::resize(esint nrows, esint ncols, esint nhalo, esint nneighbors)
{
	DataMatrixDense::resize(nrows, ncols);
	DataDistributed::resize(nhalo, info::mpi::size, nneighbors);
}

void MatrixDenseDistributed::structureUpdated()
{

}

void MatrixDenseDistributed::swap(Matrix *other)
{
	MatrixDenseDistributed *_other = other->downcast<MatrixDenseDistributed>();
	Matrix::_swap(other);
	DataMatrixDense::swap(_other);
	DataDistributed::swap(_other);
}

void MatrixDenseDistributed::shallowCopy(const Matrix *other)
{
	const MatrixDenseDistributed *_other = other->downcast<MatrixDenseDistributed>();
	Matrix::_assign(other);
	DataMatrixDense::shallowCopy(_other);
	DataDistributed::shallowCopy(_other);
}

void MatrixDenseDistributed::shallowCopyStructure(const Matrix *other)
{
	if (dynamic_cast<const MatrixDenseDistributed*>(other)) {
		const MatrixDenseDistributed *_other = other->downcast<MatrixDenseDistributed>();
		Matrix::_assign(other);
		DataMatrixDense::deepCopy(_other);
		DataDistributed::shallowCopy(_other);
		return;
	}
	if (dynamic_cast<const MatrixCSRDistributed*>(other)) {
		const MatrixCSRDistributed *_other = other->downcast<MatrixCSRDistributed>();
		Matrix::_assign(other);
		DataMatrixDense::resize(_other->nrows, _other->ncols);
		DataDistributed::shallowCopy(_other);
		return;
	}
	eslog::internalFailure("unsupported math operation: %s->fillData(%s).\n", name(), other->name());
}

void MatrixDenseDistributed::deepCopy(const Matrix *other)
{
	const MatrixDenseDistributed *_other = other->downcast<MatrixDenseDistributed>();
	Matrix::_assign(other);
	DataMatrixDense::deepCopy(_other);
	DataDistributed::deepCopy(_other);
}

void MatrixDenseDistributed::deepCopyStructure(const Matrix *other)
{
	const MatrixDenseDistributed *_other = other->downcast<MatrixDenseDistributed>();
	Matrix::_assign(other);
	DataMatrixDense::deepCopy(_other);
	DataDistributed::deepCopy(_other);
}

void MatrixDenseDistributed::uniformCombination(const Matrix *first, const Matrix *second, int nfirst, int nsecond)
{
	const MatrixDenseDistributed *_first = first->downcast<MatrixDenseDistributed>();
	const MatrixDenseDistributed *_second = second->downcast<MatrixDenseDistributed>();
	DataMatrixDense::uniformCombination(_first, _second, nfirst, nsecond);
	DataDistributed::uniformCombination(_first, _second, nfirst, nsecond);
	structureUpdated();
}

void MatrixDenseDistributed::fill(double value)
{
	DataMatrixDense::fill(value);
}

void MatrixDenseDistributed::fillData(const Matrix *in)
{
	const MatrixDenseDistributed *_in = in->downcast<MatrixDenseDistributed>();
	DataMatrixDense::fillValues(_in->vals);
}

void MatrixDenseDistributed::fillCombinedData(const Matrix *in,  esint roffset, esint coffset, esint rsize, esint csize, esint rsum, esint csum)
{
	DataMatrixDense::fillCombinedValues(in->downcast<MatrixDenseDistributed>(), roffset, coffset, rsize, csize, rsum, csum);
}

void MatrixDenseDistributed::apply(const Vector *in, Vector *out)
{
	eslog::internalFailure("call empty function.\n");
}

void MatrixDenseDistributed::apply(const Vectors *in, Vectors *out)
{
	eslog::internalFailure("call empty function.\n");
}

void MatrixDenseDistributed::scale(double alpha)
{
	MATH::vecScale(nrows * ncols, alpha, vals);
}

void MatrixDenseDistributed::add(double alpha, const Matrix *a)
{
	const MatrixDenseDistributed *_a = a->downcast<MatrixDenseDistributed>();
	MATH::vecAdd(nrows * ncols, vals, alpha, _a->vals);
}

void MatrixDenseDistributed::sum(double alpha, const Matrix *a, double beta, const Matrix *b)
{
	const MatrixDenseDistributed *_a = a->downcast<MatrixDenseDistributed>();
	const MatrixDenseDistributed *_b = b->downcast<MatrixDenseDistributed>();
	MATH::vecSum(nrows * ncols, vals, alpha, _a->vals, beta, _b->vals);
}

void MatrixDenseDistributed::addToCombination(double scale, const Matrix *in,  esint roffset, esint coffset, esint rsize, esint csize, esint rsum, esint csum)
{
	DataMatrixDense::addToCombination(scale, in->downcast<MatrixDenseDistributed>(), roffset, coffset, rsize, csize, rsum, csum);
}

void MatrixDenseDistributed::fillDiagonal(Vector *diagonal) const
{
	VectorDenseDistributed *_diagonal = diagonal->downcast<VectorDenseDistributed>();
	DataMatrixDense::fillDiagonal(_diagonal);
	_diagonal->scatterToUpper();
}

double MatrixDenseDistributed::norm()
{
	return MATH::vecNorm(nrows * ncols - nhalo * ncols, vals + nhalo * ncols);
}










