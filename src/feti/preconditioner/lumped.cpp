
#include "lumped.h"

#include "esinfo/eslog.hpp"

namespace espreso {

template <typename T>
static void _info(Lumped<T> *dual)
{
	eslog::info(" = LUMPED PRECONDITIONER PROPERTIES                                                          = \n");
	if (dual->feti->configuration.exhaustive_info) {

	}
	eslog::info(" = ----------------------------------------------------------------------------------------- = \n");
}

template <typename T>
static void _set(Lumped<T> *dual)
{

}

template <typename T>
static void _update(Lumped<T> *dual)
{

}


template <> Lumped<double>::Lumped(FETI<double> *feti): Preconditioner(feti) { _set<double>(this); }
template <> Lumped<std::complex<double> >::Lumped(FETI<std::complex<double> > *feti): Preconditioner(feti) { _set<std::complex<double> >(this); }

template <> void Lumped<double>::info() { _info<double>(this); }
template <> void Lumped<std::complex<double> >::info() { _info<std::complex<double> >(this); }

template <> void Lumped<double>::update() { _update<double>(this); }
template <> void Lumped<std::complex<double> >::update() { _update<std::complex<double> >(this); }

}
