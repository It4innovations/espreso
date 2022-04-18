
#ifndef SRC_WRAPPERS_MATH_MATRIXDENSEFETI_H_
#define SRC_WRAPPERS_MATH_MATRIXDENSEFETI_H_

#include "matrix.dense.h"
#include "matrix.feti.h"

namespace espreso {

class MatrixDenseFETI: public MatrixFETI
{
public:
	MatrixDenseFETI();
	MatrixDenseFETI(const MatrixDenseFETI &other);
	MatrixDenseFETI& operator=(const MatrixDenseFETI &other);

	MatrixDense& operator[](esint domain) { return *reinterpret_cast<MatrixDense*>(matrices[domain]); }
	const MatrixDense& operator[](esint domain) const { return *reinterpret_cast<MatrixDense*>(matrices[domain]); }

	MatrixDense* at(esint domain) { return reinterpret_cast<MatrixDense*>(matrices[domain]); }
	const MatrixDense* at(esint domain) const { return reinterpret_cast<MatrixDense*>(matrices[domain]); }

	MatrixDenseFETI* copy();

	const char* name() const { return "MatrixDenseFETI"; }

protected:
	Matrix* create();
};

}




#endif /* SRC_WRAPPERS_MATH_MATRIXDENSEFETI_H_ */
