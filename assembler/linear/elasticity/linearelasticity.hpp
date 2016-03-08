
#include "linearelasticity.h"

namespace assembler {

template <class TInput>
void LinearElasticity<TInput>::inertia(std::vector<double> &inertia)
{
	inertia.resize(3, 0);
	inertia[2] = 9.8066 * 7850; //9810.0 * 7.85e-9;
}

template <class TInput>
void LinearElasticity<TInput>::C(DenseMatrix &C, eslocal material)
{
	double ex = 2e11; //!material ? 2.1e5 : 2.1e5 * esconfig::mesh::materialDifference;
	double mi = 0.3;
	double E = ex / ((1 + mi) * (1 - 2 * mi));
	C.resize(6, 6);

	C(0, 1) = C(0, 2) = C(1, 0) = C(1, 2) = C(2, 0) = C(2, 1) = E * mi;
	C(0, 0) = C(1, 1) = C(2, 2) = E * (1.0 - mi);
	C(3, 3) = C(4, 4) = C(5, 5) = E * (0.5 - mi);
}

}
