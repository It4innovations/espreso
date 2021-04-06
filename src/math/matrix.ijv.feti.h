
#ifndef SRC_WRAPPERS_MATH_MATRIXIJVFETI_H_
#define SRC_WRAPPERS_MATH_MATRIXIJVFETI_H_

#include "matrix.ijv.h"
#include "matrix.feti.h"

namespace espreso {

class MatrixIJVFETI: public MatrixFETI
{
public:
	MatrixIJVFETI();
	MatrixIJVFETI(const MatrixIJVFETI &other);
	MatrixIJVFETI& operator=(const MatrixIJVFETI &other);

	MatrixIJV& operator[](esint domain) { return *reinterpret_cast<MatrixIJV*>(matrices[domain]); }
	const MatrixIJV& operator[](esint domain) const { return *reinterpret_cast<MatrixIJV*>(matrices[domain]); }

	MatrixIJV* at(esint domain) { return reinterpret_cast<MatrixIJV*>(matrices[domain]); }
	const MatrixIJV* at(esint domain) const { return reinterpret_cast<MatrixIJV*>(matrices[domain]); }

	MatrixIJVFETI* copy();

	const char* name() const { return "MatrixIJVFETI"; }

protected:
	Matrix* create();
};

}

#endif /* SRC_WRAPPERS_MATH_MATRIXIJVFETI_H_ */
