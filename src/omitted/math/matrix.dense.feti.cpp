
#include "matrix.dense.feti.h"
#include "vector.dense.feti.h"

using namespace espreso;

MatrixDenseFETI::MatrixDenseFETI()
{

}

MatrixDenseFETI::MatrixDenseFETI(const MatrixDenseFETI &other)
{
	deepCopy(&other);
}

MatrixDenseFETI& MatrixDenseFETI::operator=(const MatrixDenseFETI &other)
{
	if (this != &other) {
		deepCopy(&other);
	}
	return *this;
}

MatrixDenseFETI* MatrixDenseFETI::copy()
{
	return new MatrixDenseFETI();
}

Matrix* MatrixDenseFETI::create()
{
	return new MatrixDense();
}





