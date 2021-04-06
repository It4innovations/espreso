
#include "vector.h"
#include "math.h"
#include "esinfo/eslog.hpp"

#include <cstddef>

using namespace espreso;

Vector::~Vector()
{

}

Vector* Vector::shallowCopy()
{
	Vector *copy = this->copy();
	copy->shallowCopy(this);
	return copy;
}

Vector* Vector::shallowCopyStructure()
{
	Vector *copy = this->copy();
	copy->shallowCopyStructure(this);
	return copy;
}

Vector* Vector::deepCopy()
{
	Vector *copy = this->copy();
	copy->deepCopy(this);
	return copy;
}

Vector* Vector::deepCopyStructure()
{
	Vector *copy = this->copy();
	copy->deepCopyStructure(this);
	return copy;
}

void Vector::downcastFailed(const Vector *v, const Vector *targer) const
{
	eslog::error("ESPRESO internal error: cannot downcast %s into %s.\n", v->name(), targer->name());
}


Vectors::Vectors()
: nvectors(0), _holder(NULL), _vectors(NULL)
{

}

Vectors::~Vectors()
{
	if (_holder) { delete _holder; }
	for (esint n = 0; n < nvectors; ++n) {
		delete at(n);
	}
	if (_vectors) { delete[] _vectors; }
}

Vectors* Vectors::shallowCopy()
{
	Vectors *copy = this->copy();
	copy->shallowCopy(this);
	return copy;
}

Vectors* Vectors::shallowCopyStructure()
{
	Vectors *copy = this->copy();
	copy->shallowCopyStructure(this);
	return copy;
}

Vectors* Vectors::deepCopy()
{
	Vectors *copy = this->copy();
	copy->deepCopy(this);
	return copy;
}

Vectors* Vectors::deepCopyStructure()
{
	Vectors *copy = this->copy();
	copy->deepCopyStructure(this);
	return copy;
}

void Vectors::initVectors(esint nvectors)
{
	for (esint n = 0; n < this->nvectors; ++n) {
		delete at(n);
	}
	if (_vectors) { delete[] _vectors; }
	if (_holder) { delete _holder; }

	this->nvectors = nvectors;
	_holder = create();
	_vectors = new Vector*[nvectors];
	for (esint n = 0; n < nvectors; ++n) {
		_vectors[n] = create();
	}
}

void Vectors::swap(Vectors *other)
{
	_holder->swap(other->_holder);
	for (esint n = 0; n < nvectors; ++n) {
		at(n)->swap(other->at(n));
	}
}

void Vectors::shallowCopy(const Vectors *other)
{
	initVectors(other->nvectors);
	_holder->shallowCopy(other->_holder);
	for (esint n = 0; n < nvectors; ++n) {
		at(n)->shallowCopyFromHolder(_holder, n, nvectors);
	}
}

void Vectors::shallowCopyStructure(const Vectors *other)
{
	initVectors(other->nvectors);
	_holder->shallowCopyStructure(other->_holder);
	for (esint n = 0; n < nvectors; ++n) {
		at(n)->shallowCopyFromHolder(_holder, n, nvectors);
	}
}

void Vectors::deepCopy(const Vectors *other)
{
	initVectors(other->nvectors);
	_holder->deepCopy(other->_holder);
	for (esint n = 0; n < nvectors; ++n) {
		at(n)->shallowCopyFromHolder(_holder, n, nvectors);
	}
}

void Vectors::deepCopyStructure(const Vectors *other)
{
	initVectors(other->nvectors);
	_holder->deepCopyStructure(other->_holder);
	for (esint n = 0; n < nvectors; ++n) {
		at(n)->shallowCopyFromHolder(_holder, n, nvectors);
	}
}

void Vectors::uniformCombination(const Vector *first, const Vector *second, int nfirst, int nsecond)
{
	initVectors(1);
	_holder->uniformCombination(first, second, nfirst, nsecond);
	for (esint n = 0; n < nvectors; ++n) {
		at(n)->shallowCopyFromHolder(_holder, n, nvectors);
	}
}

void Vectors::uniformCombination(const Vectors *first, const Vectors *second, int nfirst, int nsecond)
{
	initVectors(first->nvectors);;
	_holder->uniformCombination(first->_holder, second->_holder, nfirst, nsecond);
	for (esint n = 0; n < nvectors; ++n) {
		at(n)->shallowCopyFromHolder(_holder, n, nvectors);
	}
}

void Vectors::fill(double value)
{
	_holder->fill(value);
}

void Vectors::fillData(const Vectors *in)
{
	for (esint n = 0; n < nvectors; ++n) {
		at(n)->fillData(in->at(n));
	}
}

void Vectors::fillCombinedValues(const Vectors *in, esint offset, esint nsize, esint sumsize)
{
	for (esint n = 0; n < nvectors; ++n) {
		at(n)->fillCombinedValues(in->at(n), offset, nsize, sumsize);
	}
}

void Vectors::fillValuesFromCombination(const Vectors *in, esint offset, esint nsize, esint sumsize)
{
	for (esint n = 0; n < nvectors; ++n) {
		at(n)->fillValuesFromCombination(in->at(n), offset, nsize, sumsize);
	}
}

void Vectors::scale(double alpha)
{
	for (esint n = 0; n < nvectors; ++n) {
		at(n)->scale(alpha);
	}
}

void Vectors::add(double alpha, const Vectors *a)
{
	for (esint n = 0; n < nvectors; ++n) {
		at(n)->add(alpha, a->at(n));
	}
}

void Vectors::sum(double alpha, const Vectors *a, double beta, const Vectors *b)
{
	for (esint n = 0; n < nvectors; ++n) {
		at(n)->sum(alpha, a->at(n), beta, b->at(n));
	}
}

void Vectors::addToCombination(double alpha, const Vectors *in, esint offset, esint nsize, esint sumsize)
{
	for (esint n = 0; n < nvectors; ++n) {
		at(n)->addToCombination(alpha, in->at(n), offset, nsize, sumsize);
	}
}

void Vectors::downcastFailed(const Vectors *v, const Vectors *target) const
{
	eslog::error("ESPRESO internal error: cannot downcast %s into %s.\n", v->name(), target->name());
}


