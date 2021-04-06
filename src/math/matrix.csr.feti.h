
#ifndef SRC_WRAPPERS_MATH_MATRIXCSRFETI_H_
#define SRC_WRAPPERS_MATH_MATRIXCSRFETI_H_

#include "matrix.csr.h"
#include "matrix.feti.h"

namespace espreso {

class MatrixDenseFETI;

class MatrixCSRFETI: public MatrixFETI
{
public:
	MatrixCSR& operator[](esint domain) { return *reinterpret_cast<MatrixCSR*>(matrices[domain]); }
	const MatrixCSR& operator[](esint domain) const { return *reinterpret_cast<MatrixCSR*>(matrices[domain]); }

	MatrixCSR* at(esint domain) { return reinterpret_cast<MatrixCSR*>(matrices[domain]); }
	const MatrixCSR* at(esint domain) const { return reinterpret_cast<MatrixCSR*>(matrices[domain]); }

	MatrixCSRFETI();
	MatrixCSRFETI(const MatrixCSRFETI &other);
	MatrixCSRFETI& operator=(const MatrixCSRFETI &other);

	MatrixCSRFETI* copy();

	void multiply(MatrixCSRFETI &A, MatrixCSRFETI &B, bool transposeA = false);
	void solve(const MatrixDenseFETI &rhs, MatrixDenseFETI &solution);
	void removeLower(MatrixType type);

	const char* name() const { return "MatrixCSRFETI"; }

protected:
	Matrix* create();
};

}

#endif /* SRC_WRAPPERS_MATH_MATRIXCSRFETI_H_ */
