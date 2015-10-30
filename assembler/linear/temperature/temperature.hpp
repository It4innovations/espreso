
#include "temperature.h"

namespace assembler {

template <MatrixComposer TMatrixComposer>
void Temperature<TMatrixComposer>::inertia(std::vector<double> &inertia)
{
	inertia.resize(1, 0);
}

template <MatrixComposer TMatrixComposer>
void Temperature<TMatrixComposer>::C(DenseMatrix &C)
{
	C.resize(3, 3);
	C(0, 0) = C(1, 1) = C(2, 2) = 1;
}

}
