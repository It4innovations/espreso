
#include "matrix.ijv.feti.h"
#include "vector.dense.feti.h"

using namespace espreso;

MatrixIJVFETI::MatrixIJVFETI()
{

}

MatrixIJVFETI::MatrixIJVFETI(const MatrixIJVFETI &other)
{
	deepCopy(&other);
}

MatrixIJVFETI& MatrixIJVFETI::operator=(const MatrixIJVFETI &other)
{
	if (this != &other) {
		deepCopy(&other);
	}
	return *this;
}

MatrixIJVFETI* MatrixIJVFETI::copy()
{
	return new MatrixIJVFETI();
}

Matrix* MatrixIJVFETI::create()
{
	return new MatrixIJV();
}


