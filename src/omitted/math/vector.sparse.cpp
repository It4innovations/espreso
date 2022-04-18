
#include "vector.sparse.h"
#include "vector.dense.h"
#include "vector.dense.feti.h"
#include "vector.dense.distributed.h"
#include "basis/containers/serializededata.h"
#include "esinfo/eslog.hpp"
#include "math.h"

#include <cstddef>

using namespace espreso;

VectorSparse::VectorSparse()
{

}

VectorSparse::VectorSparse(esint size, esint nnz)
: DataVectorSparse(size, nnz)
{

}

VectorSparse::~VectorSparse()
{

}

VectorSparse* VectorSparse::copy()
{
	return new VectorSparse();
}

void VectorSparse::swap(Vector *other)
{
	DataVectorSparse::swap(other->downcast<VectorSparse>());
}

void VectorSparse::shallowCopy(const Vector *other)
{
	DataVectorSparse::shallowCopy(other->downcast<VectorSparse>());
}

void VectorSparse::shallowCopyStructure(const Vector *other)
{
	DataVectorSparse::shallowCopyStructure(other->downcast<VectorSparse>());
}

void VectorSparse::shallowCopyFromHolder(const Vector *other, esint offset, esint nvectors)
{
	const VectorSparse *_other = other->downcast<VectorSparse>();
	DataVectorSparse::shallowCopyFromHolder(_other, offset, nvectors);
}

void VectorSparse::deepCopy(const Vector *other)
{
	DataVectorSparse::deepCopy(other->downcast<VectorSparse>());
}

void VectorSparse::deepCopyStructure(const Vector *other)
{
	DataVectorSparse::deepCopyStructure(other->downcast<VectorSparse>());
}

void VectorSparse::uniformCombination(const Vector *first, const Vector *second, int nfirst, int nsecond)
{
	DataVectorSparse::uniformCombination(first->downcast<DataVectorSparse>(), second->downcast<DataVectorSparse>(), nfirst, nsecond);
}

void VectorSparse::fill(double value)
{
	DataVectorSparse::fill(value);
}

void VectorSparse::fillData(const Vector *in)
{
	if (dynamic_cast<const VectorDense*>(in)) {
		const VectorDense *_in = dynamic_cast<const VectorDense*>(in);
		DataVectorSparse::fillDenseValues(_in->vals);
		return;
	}
	if (dynamic_cast<const VectorDenseDistributed*>(in)) {
		const VectorDenseDistributed *_in = dynamic_cast<const VectorDenseDistributed*>(in);
		DataVectorSparse::fillDenseValues(_in->vals);
		return;
	}
	if (dynamic_cast<const VectorDenseFETI*>(in)) {
		const VectorDenseFETI *_in = dynamic_cast<const VectorDenseFETI*>(in);
		auto dmap = _in->dmap->cbegin();
		for (esint i = 0, prev = 0; i < nnz; prev = indices[i++]) {
			dmap += indices[i] - prev;
			for (auto di = dmap->begin(); di != dmap->end(); ++di) {
				if (_in->ismy(di->domain)) {
					vals[i] = _in->at(di->domain - _in->doffset)->vals[di->index];
					if (_in->duplications == DataDecomposition::DUPLICATION::DUPLICATE) {
						break;
					}
				}
			}
		}
		return;
	}
	if (dynamic_cast<const VectorSparse*>(in)) {
		DataVectorSparse::fillValues(dynamic_cast<const VectorSparse*>(in)->vals);
		return;
	}
	eslog::internalFailure("unsupported math operation: VectorSparse->fillData(%s)\n", in->name());
}

void VectorSparse::fillCombinedValues(const Vector *in, esint offset, esint nsize, esint sumsize)
{
	DataVectorSparse::fillCombinedValues(in->downcast<VectorSparse>(), offset, nsize, sumsize);
}

void VectorSparse::fillValuesFromCombination(const Vector *in, esint offset, esint nsize, esint sumsize)
{
	DataVectorSparse::fillValuesFromCombination(in->downcast<VectorSparse>(), offset, nsize, sumsize);
}

void VectorSparse::scale(double alpha)
{
	MATH::vecScale(nnz, alpha, vals);
}

void VectorSparse::add(double alpha, const Vector *a)
{
	if (dynamic_cast<const VectorDense*>(a)) {
		MATH::vecAddToSparse(nnz, vals, alpha, indices, dynamic_cast<const VectorDense*>(a)->vals);
		return;
	}
	if (dynamic_cast<const VectorDenseDistributed*>(a)) {
		MATH::vecAddToSparse(nnz, vals, alpha, indices, dynamic_cast<const VectorDenseDistributed*>(a)->vals);
		return;
	}
	if (dynamic_cast<const VectorSparse*>(a)) {
		MATH::vecAdd(nnz, vals, alpha, dynamic_cast<const VectorSparse*>(a)->vals);
		return;
	}
	if (dynamic_cast<const VectorDenseFETI*>(a)) {
		const VectorDenseFETI *_a = dynamic_cast<const VectorDenseFETI*>(a);
		auto dmap = _a->dmap->cbegin();
		for (esint i = 0, prev = 0; i < nnz; prev = indices[i++]) {
			dmap += indices[i] - prev;
			for (auto di = dmap->begin(); di != dmap->end(); ++di) {
				if (_a->ismy(di->domain)) {
					vals[i] += alpha * _a->at(di->domain - _a->doffset)->vals[di->index];
					break;
				}
			}
		}
		return;
	}
	eslog::internalFailure("unsupported math operation: VectorSparse->add(%s)\n", a->name());
}

void VectorSparse::sum(double alpha, const Vector *a, double beta, const Vector *b)
{
	const VectorSparse *_a = a->downcast<VectorSparse>();
	const VectorSparse *_b = b->downcast<VectorSparse>();
	MATH::vecSum(size, vals, alpha, _a->vals, beta, _b->vals);
}

void VectorSparse::addToCombination(double alpha, const Vector *in, esint offset, esint nsize, esint sumsize)
{
	eslog::internalFailure("not implemented method\n");
}

double VectorSparse::norm()
{
	eslog::internalFailure("not implemented method\n");
	return 0;
}

double VectorSparse::max()
{
	double max = nnz ? vals[0] : 0;
	for (esint i = 1; i < nnz; ++i) {
		if (max < vals[i]) {
			max = vals[i];
		}
	}
	return max;
}

double VectorSparse::absmax()
{
	double max = nnz ? std::fabs(vals[0]) : 0;
	for (esint i = 1; i < nnz; ++i) {
		if (max < std::fabs(vals[i])) {
			max = std::fabs(vals[i]);
		}
	}
	return max;
}

double VectorSparse::dot(const Vector *other)
{
	eslog::internalFailure("not implemented method\n");
	return 0;
}

void VectorSparse::toFETI(VectorFETI *other) const
{
	if (dynamic_cast<VectorDenseFETI*>(other)) {
		VectorDenseFETI *_other = dynamic_cast<VectorDenseFETI*>(other);
		auto dmap = _other->dmap->begin();
		for (esint i = 0, prev = 0; i < nnz; prev = indices[i++]) {
			dmap += indices[i] - prev;
			for (auto di = dmap->begin(); di != dmap->end(); ++di) {
				if (_other->ismy(di->domain)) {
					_other->at(di->domain - _other->doffset)->vals[di->index] = vals[i];
				}
			}
		}
	}
}

void VectorSparse::toCombinedFETI(VectorFETI *other, esint offset, esint nsize, esint sumsize) const
{
	if (dynamic_cast<VectorDenseFETI*>(other)) {
		VectorDenseFETI *_other = dynamic_cast<VectorDenseFETI*>(other);
		auto dmap = _other->dmap->begin();
		for (esint i = 0, prev = 0; i < nnz; prev = indices[i++]) {
			dmap += indices[i] - prev;
			if (offset <= indices[i] % sumsize && indices[i] % sumsize < offset + nsize) {
				for (auto di = dmap->begin(); di != dmap->end(); ++di) {
					if (_other->ismy(di->domain)) {
						_other->at(di->domain - _other->doffset)->vals[di->index] = vals[i];
					}
				}
			}
		}
	}
}

VectorsSparse::VectorsSparse()
{
	initVectors(0);
}

VectorsSparse::~VectorsSparse()
{

}

VectorsSparse* VectorsSparse::copy()
{
	return new VectorsSparse();
}

Vector* VectorsSparse::create()
{
	return new VectorSparse();
}

void VectorsSparse::resize(esint size, esint nnz)
{
	holder()->resize(nvectors * size, nvectors * nnz);
	for (esint n = 0; n < nvectors; ++n) {
		at(n)->shallowCopyFromHolder(holder(), n, nvectors);
	}
}


