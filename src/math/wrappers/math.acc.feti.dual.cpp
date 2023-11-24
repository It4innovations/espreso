
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

#endif

template struct AccFETIDualOperator<double, Matrix_CSR>;
template struct AccFETIDualOperator<std::complex<double>, Matrix_CSR>;

}
