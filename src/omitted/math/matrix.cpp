
#include "matrix.h"
#include "esinfo/eslog.hpp"

using namespace espreso;

Matrix::Matrix()
: type(MatrixType::REAL_UNSYMMETRIC)
{

}

Matrix::Matrix(const Matrix &other)
: type(other.type)
{

}

Matrix* Matrix::shallowCopy()
{
	Matrix *copy = this->copy();
	copy->shallowCopy(this);
	return copy;
}

Matrix* Matrix::shallowCopyStructure()
{
	Matrix *copy = this->copy();
	copy->shallowCopyStructure(this);
	return copy;
}

Matrix* Matrix::deepCopy()
{
	Matrix *copy = this->copy();
	copy->deepCopy(this);
	return copy;
}

Matrix* Matrix::deepCopyStructure()
{
	Matrix *copy = this->copy();
	copy->deepCopyStructure(this);
	return copy;
}

Matrix::~Matrix()
{

}

void Matrix::downcastFailed(const Matrix *m, const Matrix *target) const
{
	eslog::internalFailure("cannot downcast %s into %s.\n", m->name(), target->name());
}

void Matrix::_assign(const Matrix *other)
{
	type = other->type;
}

void Matrix::_swap(Matrix *other)
{
	MatrixType _type = type;
	type = other->type;
	other->type = _type;
}
