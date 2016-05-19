
#include "temperature.h"

namespace espreso {

template <class TInput>
void Temperature<TInput>::inertia(std::vector<double> &inertia, const Material &material)
{
	inertia.resize(1, 0);
}

template <class TInput>
void Temperature<TInput>::C(DenseMatrix &C, const Material &material)
{
	C.resize(3, 3);
	C(0, 0) = C(1, 1) = C(2, 2) = 1;
}

}
