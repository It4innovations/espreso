
#include "vector.dense.distributed.h"
#include "data.synchronization.h"
#include "vector.sparse.h"
#include "vector.dense.h"
#include "vector.dense.feti.h"

#include "basis/containers/serializededata.h"
#include "wrappers/mpi/communication.h"
#include "esinfo/mpiinfo.h"
#include "esinfo/eslog.h"
#include "math/math.h"

#include <cstddef>
#include <cmath>

using namespace espreso;

VectorDenseDistributed::VectorDenseDistributed()
{

}

VectorDenseDistributed::VectorDenseDistributed(esint size, esint nhalo, esint nneighbors)
: DataVectorDense(size), DataDistributed(nhalo, info::mpi::size, nneighbors)
{

}

VectorDenseDistributed::~VectorDenseDistributed()
{

}

VectorDenseDistributed* VectorDenseDistributed::copy()
{
	return new VectorDenseDistributed();
}

void VectorDenseDistributed::structeUpdated()
{
	sync->init(this);
}

void VectorDenseDistributed::swap(Vector *other)
{
	VectorDenseDistributed* _other = other->downcast<VectorDenseDistributed>();
	DataVectorDense::swap(_other);
	DataDistributed::swap(_other);
}

void VectorDenseDistributed::shallowCopy(const Vector *other)
{
	const VectorDenseDistributed* _other = other->downcast<VectorDenseDistributed>();
	DataVectorDense::shallowCopyStructure(_other);
	DataDistributed::shallowCopy(_other);
}

void VectorDenseDistributed::shallowCopyStructure(const Vector *other)
{
	const VectorDenseDistributed* _other = other->downcast<VectorDenseDistributed>();
	DataVectorDense::shallowCopyStructure(_other);
	DataDistributed::shallowCopy(_other);
}

void VectorDenseDistributed::shallowCopyFromHolder(const Vector *other, esint offset, esint nvectors)
{
	const VectorDenseDistributed *_other = other->downcast<VectorDenseDistributed>();
	DataVectorDense::shallowCopyFromHolder(_other, offset, nvectors);
	DataDistributed::shallowCopy(_other);
}

void VectorDenseDistributed::deepCopy(const Vector *other)
{
	const VectorDenseDistributed* _other = other->downcast<VectorDenseDistributed>();
	DataVectorDense::deepCopy(_other);
	DataDistributed::deepCopy(_other);
}

void VectorDenseDistributed::deepCopyStructure(const Vector *other)
{
	const VectorDenseDistributed* _other = other->downcast<VectorDenseDistributed>();
	DataVectorDense::deepCopyStructure(_other);
	DataDistributed::deepCopy(_other);
}

void VectorDenseDistributed::uniformCombination(const Vector *first, const Vector *second, int nfirst, int nsecond)
{
	const VectorDenseDistributed* _first = first->downcast<VectorDenseDistributed>();
	const VectorDenseDistributed* _second = second->downcast<VectorDenseDistributed>();
	DataVectorDense::uniformCombination(_first, _second);
	DataDistributed::uniformCombination(_first, _second, nfirst, nsecond);
}

void VectorDenseDistributed::fill(double value)
{
	DataVectorDense::fill(value);
}

void VectorDenseDistributed::fillData(const Vector *in)
{
	if (dynamic_cast<const VectorDenseDistributed*>(in)) {
		DataVectorDense::fillValues(dynamic_cast<const VectorDenseDistributed*>(in)->vals);
		return;
	}

	if (dynamic_cast<const VectorDense*>(in)) {
		DataVectorDense::fillValues(dynamic_cast<const VectorDense*>(in)->vals);
		return;
	}

	if (dynamic_cast<const VectorSparse*>(in)) {
		const VectorSparse *_in = dynamic_cast<const VectorSparse*>(in);
		DataVectorDense::fillSparseValues(_in->nnz, _in->indices, _in->vals);
		return;
	}

	if (dynamic_cast<const VectorDenseFETI*>(in)) {
		const VectorDenseFETI *_in = dynamic_cast<const VectorDenseFETI*>(in);
		if (_in->duplications == DataDecomposition::DUPLICATION::DUPLICATE) {
			esint i = 0;
			for (auto dmap = _in->dmap->begin(); dmap != _in->dmap->end(); ++dmap, ++i) {
				for (auto di = dmap->begin(); di != dmap->end(); ++di) {
					if (_in->ismy(di->domain)) {
						vals[i] = _in->at(di->domain - _in->doffset)->vals[di->index];
						break;
					}
				}
			}
			return;
		}
	}
	eslog::internalFailure("unsupported math operation.\n");
}

void VectorDenseDistributed::fillCombinedValues(const Vector *in, esint offset, esint nsize, esint sumsize)
{
	if (dynamic_cast<const VectorDense*>(in)) {
		const VectorDense* _in = in->downcast<VectorDense>();
		DataVectorDense::fillCombinedValues(_in, offset, nsize, sumsize);
		return;
	}
	if (dynamic_cast<const VectorDenseDistributed*>(in)) {
		const VectorDenseDistributed* _in = in->downcast<VectorDenseDistributed>();
		DataVectorDense::fillCombinedValues(_in, offset, nsize, sumsize);
		return;
	}

	eslog::internalFailure("unsupported math operation.\n");
}

void VectorDenseDistributed::fillValuesFromCombination(const Vector *in, esint offset, esint nsize, esint sumsize)
{
	const VectorDenseDistributed* _in = in->downcast<VectorDenseDistributed>();
	DataVectorDense::fillValuesFromCombination(_in, offset, nsize, sumsize);
}

void VectorDenseDistributed::scale(double alpha)
{
	MATH::vecScale(size, alpha, vals);
}

void VectorDenseDistributed::add(double alpha, const Vector *a)
{
	const VectorDenseDistributed* _a = a->downcast<VectorDenseDistributed>();
	MATH::vecAdd(size, vals, alpha, _a->vals);
}

void VectorDenseDistributed::sum(double alpha, const Vector *a, double beta, const Vector *b)
{
	const VectorDenseDistributed* _a = a->downcast<VectorDenseDistributed>();
	const VectorDenseDistributed* _b = b->downcast<VectorDenseDistributed>();
	MATH::vecSum(size, vals, alpha, _a->vals, beta, _b->vals);
}

void VectorDenseDistributed::addToCombination(double alpha, const Vector *in, esint offset, esint nsize, esint sumsize)
{
	const VectorDenseDistributed* _in = in->downcast<VectorDenseDistributed>();
	DataVectorDense::addToCombination(alpha, _in, offset, nsize, sumsize);
}

double VectorDenseDistributed::norm()
{
	double res, norm = MATH::vecDot(size - nhalo, vals + nhalo);
	Communication::allReduce(&norm, &res, 1, MPI_DOUBLE, MPI_SUM);
	return std::sqrt(res);
}

double VectorDenseDistributed::max()
{
	double max = size ? vals[nhalo] : 0;
	for (esint i = nhalo + 1; i < size; ++i) {
		if (max < vals[i]) {
			max = vals[i];
		}
	}
	double gmax;
	Communication::allReduce(&max, &gmax, 1, MPI_DOUBLE, MPI_MAX);
	return gmax;
}

double VectorDenseDistributed::absmax()
{
	double max = size ? std::fabs(vals[nhalo]) : 0;
	for (esint i = nhalo + 1; i < size; ++i) {
		if (max < std::fabs(vals[i])) {
			max = std::fabs(vals[i]);
		}
	}
	double gmax;
	Communication::allReduce(&max, &gmax, 1, MPI_DOUBLE, MPI_MAX);
	return gmax;
}

double VectorDenseDistributed::dot(const Vector *other)
{
	double res, norm = MATH::vecDot(size - nhalo, vals + nhalo, other->downcast<VectorDenseDistributed>()->vals + nhalo);
	Communication::allReduce(&norm, &res, 1, MPI_DOUBLE, MPI_SUM);
	return res;
}

void VectorDenseDistributed::gatherFromUpper()
{
	sync->gatherFromUpper(this);
}

void VectorDenseDistributed::scatterToUpper()
{
	sync->scatterToUpper(this);
}

VectorsDenseDistributed::VectorsDenseDistributed()
{
	initVectors(0);
}

VectorsDenseDistributed::~VectorsDenseDistributed()
{

}

VectorsDenseDistributed* VectorsDenseDistributed::copy()
{
	return new VectorsDenseDistributed();
}

Vector* VectorsDenseDistributed::create()
{
	return new VectorDenseDistributed();
}

void VectorsDenseDistributed::resize(esint size, esint nhalo, esint nneighbors)
{
	holder()->DataVectorDense::resize(nvectors * size);
	holder()->DataDistributed::resize(nhalo, info::mpi::size, nneighbors);
	for (esint n = 0; n < nvectors; ++n) {
		at(n)->shallowCopyFromHolder(holder(), n, nvectors);
	}
}

void VectorsDenseDistributed::structureUpdated()
{
	holder()->structeUpdated();
}

void VectorsDenseDistributed::fillDistribution(esint *halo, esint *distribution, int *neighbors)
{
	holder()->fillDistribution(halo, distribution, neighbors);
}

void VectorsDenseDistributed::gatherFromUpper()
{
	for (esint n = 0; n < nvectors; ++n) {
		at(n)->gatherFromUpper();
	}
}

void VectorsDenseDistributed::scatterToUpper()
{
	for (esint n = 0; n < nvectors; ++n) {
		at(n)->scatterToUpper();
	}
}


