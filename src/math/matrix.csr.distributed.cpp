
#include "matrix.csr.distributed.h"
#include "matrix.csr.h"
#include "vector.dense.h"
#include "vector.dense.distributed.h"
#include "data.synchronization.h"
#include "data.mv.h"
#include "math.h"

#include "esinfo/eslog.h"
#include "esinfo/mpiinfo.h"
#include "basis/utilities/communication.h"

#include <cstddef>

using namespace espreso;

MatrixCSRDistributed::MatrixCSRDistributed()
{

}

MatrixCSRDistributed::MatrixCSRDistributed(esint nrows, esint ncols, esint nnz, esint nhalo, esint nneighbors)
: DataMatrixCSR(nrows, ncols, nnz), DataDistributed(nhalo, info::mpi::size, nneighbors)
{
	structureUpdated();
}

MatrixCSRDistributed::MatrixCSRDistributed(const MatrixCSRDistributed &other)
{
	deepCopy(&other);
}

MatrixCSRDistributed& MatrixCSRDistributed::operator=(const MatrixCSRDistributed &other)
{
	if (this != &other) {
		deepCopy(&other);
	}
	return *this;
}

MatrixCSRDistributed::~MatrixCSRDistributed()
{

}


MatrixCSRDistributed* MatrixCSRDistributed::copy()
{
	return new MatrixCSRDistributed();;
}

void MatrixCSRDistributed::resize(esint nrows, esint ncols, esint nnz, esint nhalo, esint nneighbors)
{
	DataMatrixCSR::resize(nrows, ncols, nnz);
	DataDistributed::resize(nhalo, info::mpi::size, nneighbors);
}

void MatrixCSRDistributed::structureUpdated()
{
	mv->init(this);
	sync->init(this);
}

void MatrixCSRDistributed::swap(Matrix *other)
{
	MatrixCSRDistributed *_other = other->downcast<MatrixCSRDistributed>();
	Matrix::_swap(other);
	DataMatrixCSR::swap(_other);
	DataDistributed::swap(_other);
}

void MatrixCSRDistributed::shallowCopy(const Matrix *other)
{
	const MatrixCSRDistributed *_other = other->downcast<MatrixCSRDistributed>();
	Matrix::_assign(other);
	DataMatrixCSR::shallowCopy(_other);
	DataDistributed::shallowCopy(_other);
}

void MatrixCSRDistributed::shallowCopyStructure(const Matrix *other)
{
	const MatrixCSRDistributed *_other = other->downcast<MatrixCSRDistributed>();
	Matrix::_assign(other);
	DataMatrixCSR::shallowCopyStructure(_other);
	DataDistributed::shallowCopy(_other);
}

void MatrixCSRDistributed::deepCopy(const Matrix *other)
{
	const MatrixCSRDistributed *_other = other->downcast<MatrixCSRDistributed>();
	Matrix::_assign(other);
	DataMatrixCSR::deepCopy(_other);
	DataDistributed::deepCopy(_other);
}

void MatrixCSRDistributed::deepCopyStructure(const Matrix *other)
{
	const MatrixCSRDistributed *_other = other->downcast<MatrixCSRDistributed>();
	Matrix::_assign(other);
	DataMatrixCSR::deepCopyStructure(_other);
	DataDistributed::deepCopy(_other);
}

void MatrixCSRDistributed::uniformCombination(const Matrix *first, const Matrix *second, int nfirst, int nsecond)
{
	const MatrixCSRDistributed *_first = first->downcast<MatrixCSRDistributed>();
	const MatrixCSRDistributed *_second = second->downcast<MatrixCSRDistributed>();
	DataMatrixCSR::uniformCombination(_first, _second, nfirst, nsecond);
	DataDistributed::uniformCombination(_first, _second, nfirst, nsecond);
	structureUpdated();
}

void MatrixCSRDistributed::fill(double value)
{
	DataMatrixCSR::fill(value);
}

void MatrixCSRDistributed::fillData(const Matrix *in)
{
	const MatrixCSRDistributed *_in = in->downcast<MatrixCSRDistributed>();
	DataMatrixCSR::fillValues(_in->nnz, _in->vals);
}

void MatrixCSRDistributed::fillCombinedData(const Matrix *in, esint roffset, esint coffset, esint nsize, esint sumsize)
{
	DataMatrixCSR::fillCombinedValues(in->downcast<MatrixCSRDistributed>(), roffset, coffset, nsize, sumsize);
}

void MatrixCSRDistributed::apply(const Vector *in, Vector *out)
{
	const VectorDenseDistributed *_in = in->downcast<VectorDenseDistributed>();
	VectorDenseDistributed *_out = out->downcast<VectorDenseDistributed>();
	mv->apply(this, _in, _out);
}

void MatrixCSRDistributed::apply(const Vectors *in, Vectors *out)
{
	for (esint n = 0; n < in->nvectors; ++n) {
		apply(in->at(n), out->at(n));
	}
}

void MatrixCSRDistributed::scale(double alpha)
{
	MATH::vecScale(nnz, alpha, vals);
}

void MatrixCSRDistributed::add(double alpha, const Matrix *a)
{
	const MatrixCSRDistributed *_a = a->downcast<MatrixCSRDistributed>();
	MATH::vecAdd(nnz, vals, alpha, _a->vals);
}

void MatrixCSRDistributed::sum(double alpha, const Matrix *a, double beta, const Matrix *b)
{
	const MatrixCSRDistributed *_a = a->downcast<MatrixCSRDistributed>();
	const MatrixCSRDistributed *_b = b->downcast<MatrixCSRDistributed>();
	MATH::vecSum(nnz, vals, alpha, _a->vals, beta, _b->vals);
}

void MatrixCSRDistributed::addToCombination(double scale, const Matrix *in, esint roffset, esint coffset, esint nsize, esint sumsize)
{
	DataMatrixCSR::addToCombination(scale, in->downcast<MatrixCSRDistributed>(), roffset, coffset, nsize, sumsize);
}

void MatrixCSRDistributed::fillDiagonal(Vector *diagonal) const
{
	VectorDenseDistributed *_diagonal = diagonal->downcast<VectorDenseDistributed>();
	DataMatrixCSR::fillDiagonal(_diagonal);
	_diagonal->scatterToUpper();
}

double MatrixCSRDistributed::norm()
{
	return MATH::vecNorm(nnz - nhalo, vals + nhalo);
}

void MatrixCSRDistributed::gatherFromUpper()
{
	sync->gatherFromUpper(this);
}










