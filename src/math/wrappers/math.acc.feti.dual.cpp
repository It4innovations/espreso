
#include "math.acc.feti.dual.h"
#include "esinfo/eslog.h"

#include <complex>

namespace espreso {

#ifndef HAVE_ROCM

template <template <typename, typename, typename> class Matrix, typename T, typename I>
AccFETIDualOperator<Matrix, T, I>::AccFETIDualOperator(int rank)
: rank(rank), _acc(nullptr)
{
    eslog::error("calling of empty AccFETIDualOperator wrapper.\n");
}


template <template <typename, typename, typename> class Matrix, typename T, typename I>
AccFETIDualOperator<Matrix, T, I>::~AccFETIDualOperator()
{

}

template <template <typename, typename, typename> class Matrix, typename T, typename I>
void AccFETIDualOperator<Matrix, T, I>::set(const std::vector<MatrixType> &K, const std::vector<MatrixType> &B)
{

}

template <template <typename, typename, typename> class Matrix, typename T, typename I>
void AccFETIDualOperator<Matrix, T, I>::update(const std::vector<MatrixType> &K)
{

}

template <template <typename, typename, typename> class Matrix, typename T, typename I>
void AccFETIDualOperator<Matrix, T, I>::apply(const Vector_Dual<T> &x, Vector_Dual<T> &y, const std::vector<std::vector<int> > & D2C)
{

}

#endif

template struct AccFETIDualOperator<Matrix_CSR, double, int>;
template struct AccFETIDualOperator<Matrix_CSR, std::complex<double>, int>;

}
