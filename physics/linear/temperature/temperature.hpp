
#include "temperature.h"

namespace physics {

template <MatrixComposer TMatrixComposer>
void Temperature<TMatrixComposer>::init()
{
	SparseSolver ss;
	this->K(ss, 0);
}

template <MatrixComposer TMatrixComposer>
void Temperature<TMatrixComposer>::pre_solve_update()
{

}

template <MatrixComposer TMatrixComposer>
void Temperature<TMatrixComposer>::post_solve_update()
{

}

template <MatrixComposer TMatrixComposer>
void Temperature<TMatrixComposer>::solve()
{

}

template <MatrixComposer TMatrixComposer>
void Temperature<TMatrixComposer>::finalize()
{

}

template <MatrixComposer TMatrixComposer>
void Temperature<TMatrixComposer>::C(DenseMatrix &C)
{
	C.resize(3, 3);
	C(0, 0) = C(1, 1) = C(2, 2) = 1;
}

}
