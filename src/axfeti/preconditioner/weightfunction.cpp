
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


template <> WeightFunction<double>::WeightFunction(AX_FETI<double> *feti): Preconditioner(feti) { _set<double>(this); }
template <> WeightFunction<std::complex<double> >::WeightFunction(AX_FETI<std::complex<double> > *feti): Preconditioner(feti) { _set<std::complex<double> >(this); }

template <> void WeightFunction<double>::info() { _info<double>(this); }
template <> void WeightFunction<std::complex<double> >::info() { _info<std::complex<double> >(this); }

template <> void WeightFunction<double>::update() { _update<double>(this); }
template <> void WeightFunction<std::complex<double> >::update() { _update<std::complex<double> >(this); }

}
