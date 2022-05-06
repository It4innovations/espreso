
#include "weightfunction.h"

#include "esinfo/eslog.hpp"

namespace espreso {

template <typename T>
static void _info(WeightFunction<T> *dual)
{
	eslog::info(" = WEIGHT FUNCTION PRECONDITIONING PROPERTIES                                                = \n");
	if (dual->feti->configuration.exhaustive_info) {

	}
	eslog::info(" = ----------------------------------------------------------------------------------------- = \n");
}

template <typename T>
static void _set(WeightFunction<T> *dual)
{

}

template <typename T>
static void _update(WeightFunction<T> *dual)
{

}

template <typename T>
static void _apply(WeightFunction<T> *dual)
{

}

template <> WeightFunction<double>::WeightFunction(FETI<double> *feti): Preconditioner(feti) { _set<double>(this); }
template <> WeightFunction<std::complex<double> >::WeightFunction(FETI<std::complex<double> > *feti): Preconditioner(feti) { _set<std::complex<double> >(this); }

template <> void WeightFunction<double>::info() { _info<double>(this); }
template <> void WeightFunction<std::complex<double> >::info() { _info<std::complex<double> >(this); }

template <> void WeightFunction<double>::update() { _update<double>(this); }
template <> void WeightFunction<std::complex<double> >::update() { _update<std::complex<double> >(this); }

template <> void WeightFunction<double>::apply(const Vector_Dual<double> &x, Vector_Dual<double> &y) { _apply(this); }
template <> void WeightFunction<std::complex<double> >::apply(const Vector_Dual<std::complex<double> > &x, Vector_Dual<std::complex<double> > &y) { _apply(this); }

}
