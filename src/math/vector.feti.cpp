
#include "vector.feti.h"
#include "vector.dense.feti.h"
#include "vector.sparse.feti.h"
#include "esinfo/eslog.h"

#include <cstddef>

using namespace espreso;

VectorFETI::VectorFETI()
: DataDecomposition(DataDecomposition::DUPLICATION::DUPLICATE), domains(0), vectors(NULL)
{

}

VectorFETI::~VectorFETI()
{
	for (esint d = 0; d < domains; ++d) {
		delete at(d);
	}
	delete[] vectors;
}

VectorFETI* VectorFETI::copy()
{
	eslog::internalFailure("VectorFETI cannot be used.\n");
	return NULL;
}

Vector* VectorFETI::create()
{
	eslog::internalFailure("cannot create instance of VectorFETI.\n");
	return NULL;
}

void VectorFETI::setDuplications(DataDecomposition::DUPLICATION duplications)
{
	this->duplications = duplications;
}

void VectorFETI::initDomains(DataDecomposition::DUPLICATION duplications, esint domains)
{
	this->duplications = duplications;
	for (esint d = 0; d < this->domains; ++d) {
		delete at(d);
	}
	if (vectors) { delete[] vectors; }

	this->domains = domains;
	vectors = new Vector*[domains];
	#pragma omp parallel for
	for (esint d = 0; d < domains; ++d) {
		vectors[d] = create();
	}
}

void VectorFETI::swap(Vector *other)
{
	VectorFETI *_other = other->downcast<VectorFETI>();
	#pragma omp parallel for
	for (esint d = 0; d < domains; ++d) {
		at(d)->swap(_other->at(d));
	}
	DataDecomposition::swap(_other);
}

void VectorFETI::shallowCopy(const Vector *other)
{
	const VectorFETI *_other = other->downcast<VectorFETI>();
	initDomains(_other->duplications, _other->domains);
	#pragma omp parallel for
	for (esint d = 0; d < domains; ++d) {
		at(d)->shallowCopy(_other->at(d));
	}
	DataDecomposition::shallowCopy(_other);
}

void VectorFETI::shallowCopyStructure(const Vector *other)
{
	const VectorFETI *_other = other->downcast<VectorFETI>();
	initDomains(_other->duplications, _other->domains);
	#pragma omp parallel for
	for (esint d = 0; d < domains; ++d) {
		at(d)->shallowCopyStructure(_other->at(d));
	}
	DataDecomposition::shallowCopyStructure(_other);
}

void VectorFETI::shallowCopyFromHolder(const Vector *other, esint offset, esint nvectors)
{
	const VectorFETI *_other = other->downcast<VectorFETI>();
	initDomains(_other->duplications, _other->domains);
	#pragma omp parallel for
	for (esint d = 0; d < domains; ++d) {
		at(d)->shallowCopyFromHolder(_other->at(d), offset, nvectors);
	}
	DataDecomposition::shallowCopy(_other);
}

void VectorFETI::deepCopy(const Vector *other)
{
	const VectorFETI *_other = other->downcast<VectorFETI>();
	initDomains(_other->duplications, _other->domains);
	#pragma omp parallel for
	for (esint d = 0; d < domains; ++d) {
		at(d)->deepCopy(_other->at(d));
	}
	DataDecomposition::deepCopy(_other);
}

void VectorFETI::deepCopyStructure(const Vector *other)
{
	const VectorFETI *_other = other->downcast<VectorFETI>();
	initDomains(_other->duplications, _other->domains);
	#pragma omp parallel for
	for (esint d = 0; d < domains; ++d) {
		at(d)->deepCopyStructure(_other->at(d));
	}
	DataDecomposition::deepCopyStructure(_other);
}

void VectorFETI::uniformCombination(const Vector *first, const Vector *second, int nfirst, int nsecond)
{
	const VectorFETI *_first = first->downcast<VectorFETI>();
	const VectorFETI *_second = second->downcast<VectorFETI>();
	initDomains(_first->duplications, _first->domains);
	#pragma omp parallel for
	for (esint d = 0; d < domains; ++d) {
		at(d)->uniformCombination(_first->at(d), _second->at(d), nfirst, nsecond);
	}
	DataDecomposition::uniformCombination(_first, _second, nfirst, nsecond);
}

void VectorFETI::fill(double value)
{
	#pragma omp parallel for
	for (esint d = 0; d < domains; ++d) {
		at(d)->fill(value);
	}
}

void VectorFETI::fillData(const Vector *in)
{
	if (dynamic_cast<const VectorFETI*>(in)) {
		if (dynamic_cast<const VectorDenseFETI*>(in) && dynamic_cast<VectorDenseFETI*>(this)) {
			const VectorDenseFETI* _in = dynamic_cast<const VectorDenseFETI*>(in);
			VectorDenseFETI* _self = dynamic_cast<VectorDenseFETI*>(this);
			bool same = _in->domains == _self->domains;
			for (esint d = 0; same && d < _self->domains; ++d) {
				same &= _in->at(d)->size == _self->at(d)->size;
			}
			if (same) {
				#pragma omp parallel for
				for (esint d = 0; d < domains; ++d) {
					at(d)->fillData(_in->at(d));
				}
			} else {
				dynamic_cast<const VectorDenseFETI*>(in)->toFETI(_self);
			}
			return;
		}
		return;
	}
	if (dynamic_cast<const VectorDense*>(in)) {
		dynamic_cast<const VectorDense*>(in)->toFETI(this);
		return;
	}
	if (dynamic_cast<const VectorSparse*>(in)) {
		dynamic_cast<const VectorSparse*>(in)->toFETI(this);
		return;
	}
	eslog::internalFailure("unsupported math operation.\n");
}

void VectorFETI::fillCombinedValues(const Vector *in, esint offset, esint nsize, esint sumsize)
{
	if (dynamic_cast<const VectorFETI*>(in)) {
		const VectorFETI *_in = in->downcast<VectorFETI>();
		#pragma omp parallel for
		for (esint d = 0; d < domains; ++d) {
			at(d)->fillCombinedValues(_in->at(d), offset, nsize, sumsize);
		}
		return;
	}
	if (dynamic_cast<const VectorDense*>(in)) {
		dynamic_cast<const VectorDense*>(in)->toCombinedFETI(this, offset, nsize, sumsize);
		return;
	}
	if (dynamic_cast<const VectorSparse*>(in)) {
		dynamic_cast<const VectorSparse*>(in)->toCombinedFETI(this, offset, nsize, sumsize);
		return;
	}
	eslog::internalFailure("unsupported math operation.\n");
}

void VectorFETI::fillValuesFromCombination(const Vector *in, esint offset, esint nsize, esint sumsize)
{
	const VectorFETI *_in = in->downcast<VectorFETI>();
	#pragma omp parallel for
	for (esint d = 0; d < domains; ++d) {
		at(d)->fillValuesFromCombination(_in->at(d), offset, nsize, sumsize);
	}
}

void VectorFETI::scale(double alpha)
{
	#pragma omp parallel for
	for (esint d = 0; d < domains; ++d) {
		at(d)->scale(alpha);
	}
}

void VectorFETI::add(double alpha, const Vector *a)
{
	const VectorFETI *_a = a->downcast<VectorFETI>();
	#pragma omp parallel for
	for (esint d = 0; d < domains; ++d) {
		at(d)->add(alpha, _a->at(d));
	}
}

void VectorFETI::sum(double alpha, const Vector *a, double beta, const Vector *b)
{
	const VectorFETI *_a = a->downcast<VectorFETI>();
	const VectorFETI *_b = b->downcast<VectorFETI>();
	#pragma omp parallel for
	for (esint d = 0; d < domains; ++d) {
		at(d)->sum(alpha, _a->at(d), beta, _b->at(d));
	}
}

void VectorFETI::addToCombination(double alpha, const Vector *in, esint offset, esint nsize, esint sumsize)
{
	const VectorFETI *_in = in->downcast<VectorFETI>();
	#pragma omp parallel for
	for (esint d = 0; d < domains; ++d) {
		at(d)->addToCombination(alpha, _in->at(d), offset, nsize, sumsize);
	}
}

double VectorFETI::norm()
{
	eslog::internalFailure("call empty function.\n");
	return 0;
}

double VectorFETI::max()
{
	eslog::internalFailure("call empty function.\n");
	return 0;
}

double VectorFETI::absmax()
{
	eslog::internalFailure("call empty function.\n");
	return 0;
}

double VectorFETI::dot(const Vector *other)
{
	eslog::internalFailure("call empty function.\n");
	return 0;
}

void VectorFETI::fromFETI(VectorFETI *other) const
{
	eslog::internalFailure("call empty function.\n");
}

VectorsFETI::VectorsFETI()
{

}

VectorsFETI::~VectorsFETI()
{

}

VectorsFETI* VectorsFETI::copy()
{
	eslog::internalFailure("cannot copy VectorsFETI.\n");
	return NULL;
}

void VectorsFETI::setDuplications(DataDecomposition::DUPLICATION duplications)
{
	_holder->downcast<VectorFETI>()->setDuplications(duplications);
	for (esint n = 0; n < nvectors; ++n) {
		at(n)->downcast<VectorFETI>()->setDuplications(duplications);
	}
}

Vector* VectorsFETI::create()
{
	eslog::internalFailure("cannot create VectorsFETI.\n");
	return NULL;
}

