
#include "math.acc.feti.dual.h"
#include "esinfo/eslog.h"

#include <complex>

namespace espreso {

#ifndef HAVE_ROCM

template <typename T, template <typename> class Matrix>
AccFETIDualOperator<T, Matrix>::AccFETIDualOperator()
: _acc(nullptr)
{
	eslog::error("calling of empty AccFETIDualOperator wrapper.\n");
}


template <typename T, template <typename> class Matrix>
AccFETIDualOperator<T, Matrix>::~AccFETIDualOperator()
{

}

template <typename T, template <typename> class Matrix>
void AccFETIDualOperator<T, Matrix>::set(const std::vector<Matrix_CSR<T> > &K, const std::vector<Matrix_CSR<T> > &B)
{

}

template <typename T, template <typename> class Matrix>
void AccFETIDualOperator<T, Matrix>::update(const std::vector<Matrix_CSR<T> > &K)
{

}

template <typename T, template <typename> class Matrix>
void AccFETIDualOperator<T, Matrix>::apply(const Vector_Dual<T> &x, Vector_Dual<T> &y, const std::vector<std::vector<int> > & D2C)
{

}

#endif

template struct AccFETIDualOperator<double, Matrix_CSR>;
template struct AccFETIDualOperator<std::complex<double>, Matrix_CSR>;

}
