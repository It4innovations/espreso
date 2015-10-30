
#include "linearelasticity.h"

namespace assembler {

template <MatrixComposer TMatrixComposer>
void LinearElasticity<TMatrixComposer>::inertia(std::vector<double> &inertia)
{
	inertia.resize(3, 0);
	inertia[2] = 9810.0 * 7.85e-9;
}

template <MatrixComposer TMatrixComposer>
void LinearElasticity<TMatrixComposer>::C(DenseMatrix &C)
{
	std::vector<double> inertia(3, 0.0);
	inertia[2] = 9810.0 * 7.85e-9;

	double ex = 2.1e5;
	double mi = 0.3;
	double E = ex / ((1 + mi) * (1 - 2 * mi));
	C.resize(6, 6);

	C(0, 1) = C(0, 2) = C(1, 0) = C(1, 2) = C(2, 0) = C(2, 1) = E * mi;
	C(0, 0) = C(1, 1) = C(2, 2) = E * (1.0 - mi);
	C(3, 3) = C(4, 4) = C(5, 5) = E * (0.5 - mi);
}

}
