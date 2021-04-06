
#include "vector.sparse.distributed.h"

#include "basis/utilities/communication.h"
#include "esinfo/mpiinfo.h"
#include "esinfo/eslog.h"
#include "math.h"

#include <cstddef>
#include <cmath>

using namespace espreso;

VectorSparseDistributed::VectorSparseDistributed()
{

}

VectorSparseDistributed::VectorSparseDistributed(esint size, esint nnz, esint nhalo, esint nneighbors)
: DataVectorSparse(size, nnz), DataDistributed(nhalo, info::mpi::size, nneighbors)
{

}

VectorSparseDistributed::~VectorSparseDistributed()
{

}

VectorSparseDistributed* VectorSparseDistributed::copy()
{
	return new VectorSparseDistributed();
}

void VectorSparseDistributed::swap(Vector *other)
{
	VectorSparseDistributed *_other = other->downcast<VectorSparseDistributed>();
	DataVectorSparse::swap(_other);
	DataDistributed::swap(_other);
}

void VectorSparseDistributed::shallowCopy(const Vector *other)
{
	const VectorSparseDistributed *_other = other->downcast<VectorSparseDistributed>();
	DataVectorSparse::shallowCopy(_other);
	DataDistributed::shallowCopy(_other);
}

void VectorSparseDistributed::shallowCopyStructure(const Vector *other)
{
	const VectorSparseDistributed *_other = other->downcast<VectorSparseDistributed>();
	DataVectorSparse::shallowCopyStructure(_other);
	DataDistributed::shallowCopy(_other);
}

void VectorSparseDistributed::shallowCopyFromHolder(const Vector *other, esint offset, esint nvectors)
{
	const VectorSparseDistributed *_other = other->downcast<VectorSparseDistributed>();
	DataVectorSparse::shallowCopyFromHolder(_other, offset, nvectors);
	DataDistributed::shallowCopy(_other);
}

void VectorSparseDistributed::deepCopy(const Vector *other)
{
	const VectorSparseDistributed *_other = other->downcast<VectorSparseDistributed>();
	DataVectorSparse::deepCopy(_other);
	DataDistributed::deepCopy(_other);
}

void VectorSparseDistributed::deepCopyStructure(const Vector *other)
{
	const VectorSparseDistributed *_other = other->downcast<VectorSparseDistributed>();
	DataVectorSparse::deepCopyStructure(_other);
	DataDistributed::deepCopy(_other);
}

void VectorSparseDistributed::uniformCombination(const Vector *first, const Vector *second, int nfirst, int nsecond)
{
	const VectorSparseDistributed* _first = first->downcast<VectorSparseDistributed>();
	const VectorSparseDistributed* _second = second->downcast<VectorSparseDistributed>();
	DataVectorSparse::uniformCombination(_first, _second, nfirst, nsecond);
	DataDistributed::uniformCombination(_first, _second, nfirst, nsecond);
}

void VectorSparseDistributed::fill(double value)
{
	DataVectorSparse::fill(value);
}

void VectorSparseDistributed::fillData(const Vector *in)
{
	const VectorSparseDistributed* _in = in->downcast<VectorSparseDistributed>();
	DataVectorSparse::fillValues(_in->vals);
}

void VectorSparseDistributed::fillCombinedValues(const Vector *in, esint offset, esint nsize, esint sumsize)
{
	const VectorSparseDistributed* _in = in->downcast<VectorSparseDistributed>();
	DataVectorSparse::fillCombinedValues(_in, offset, nsize, sumsize);
}

void VectorSparseDistributed::fillValuesFromCombination(const Vector *in, esint offset, esint nsize, esint sumsize)
{
	const VectorSparseDistributed* _in = in->downcast<VectorSparseDistributed>();
	DataVectorSparse::fillValuesFromCombination(_in, offset, nsize, sumsize);
}

void VectorSparseDistributed::scale(double alpha)
{
	MATH::vecScale(nnz, alpha, vals);
}

void VectorSparseDistributed::add(double alpha, const Vector *a)
{
	const VectorSparseDistributed* _a = a->downcast<VectorSparseDistributed>();
	MATH::vecAdd(nnz, vals, alpha, _a->vals);
}

void VectorSparseDistributed::sum(double alpha, const Vector *a, double beta, const Vector *b)
{
	const VectorSparseDistributed* _a = a->downcast<VectorSparseDistributed>();
	const VectorSparseDistributed* _b = b->downcast<VectorSparseDistributed>();
	MATH::vecSum(nnz, vals, alpha, _a->vals, beta, _b->vals);
}

void VectorSparseDistributed::addToCombination(double alpha, const Vector *in, esint offset, esint nsize, esint sumsize)
{
	const VectorSparseDistributed* _in = in->downcast<VectorSparseDistributed>();
	DataVectorSparse::addToCombination(alpha, _in, offset, nsize, sumsize);
}

double VectorSparseDistributed::norm()
{
	double res, norm = MATH::vecDot(nnz - nhalo, vals + nhalo);
	Communication::allReduce(&norm, &res, 1, MPI_DOUBLE, MPI_SUM);
	return std::sqrt(res);
}

double VectorSparseDistributed::max()
{
	eslog::error("ESPRESO internal error: call empty function.");
	return 0;
}

double VectorSparseDistributed::absmax()
{
	eslog::error("ESPRESO internal error: call empty function.");
	return 0;
}

double VectorSparseDistributed::dot(const Vector *other)
{
	double res, norm = MATH::vecDot(nnz - nhalo, vals + nhalo, other->downcast<VectorSparseDistributed>()->vals + nhalo);
	Communication::allReduce(&norm, &res, 1, MPI_DOUBLE, MPI_SUM);
	return res;
}

void VectorSparseDistributed::averageDuplications()
{
	eslog::error("ESPRESO internal error: call empty function.\n");
//	DataVectorDistributed::synchronize(vals);
}

VectorsSparseDistributed::VectorsSparseDistributed()
{
	initVectors(0);
}

VectorsSparseDistributed::~VectorsSparseDistributed()
{

}

VectorsSparseDistributed* VectorsSparseDistributed::copy()
{
	return new VectorsSparseDistributed();
}

Vector* VectorsSparseDistributed::create()
{
	return new VectorSparseDistributed();
}

void VectorsSparseDistributed::resize(esint nvectors, esint nnz, esint nhalo, esint nneighbors)
{
	eslog::error("ESPRESO internal error: call empty function.\n");
//	holder.DataVectorSparse::resize(nvectors * nnz);
//	holder.DataVectorDistributed::resize(nhalo);
//	_update(nvectors, nnz, nhalo);
}

void VectorsSparseDistributed::fillDistribution(esint nhalo, esint *halo, esint *distribution, int *neighbors)
{
	eslog::error("ESPRESO internal error: call empty function.\n");
//	holder.fillDistribution(nhalo, halo, distribution);
//	_update(nvectors, nnz, nhalo);
}


