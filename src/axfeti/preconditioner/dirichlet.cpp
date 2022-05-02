
#include "dirichlet.h"

#include "esinfo/eslog.hpp"

using namespace espreso;

template <typename T>
static void _info(Dirichlet<T> *dual)
{
	eslog::info(" = DIRICHLET PRECONDITIONER PROPERTIES                                                       = \n");
	if (dual->feti->configuration.exhaustive_info) {

	}
	eslog::info(" = ----------------------------------------------------------------------------------------- = \n");
}

template <typename T>
static void _set(Dirichlet<T> *dual)
{

}

template <typename T>
static void _update(Dirichlet<T> *dual)
{

}


template <> Dirichlet<double>::Dirichlet(AX_FETI<double> *feti): Preconditioner(feti) { _set<double>(this); }
template <> Dirichlet<std::complex<double> >::Dirichlet(AX_FETI<std::complex<double> > *feti): Preconditioner(feti) { _set<std::complex<double> >(this); }

template <> void Dirichlet<double>::info() { _info<double>(this); }
template <> void Dirichlet<std::complex<double> >::info() { _info<std::complex<double> >(this); }

template <> void Dirichlet<double>::update() { _update<double>(this); }
template <> void Dirichlet<std::complex<double> >::update() { _update<std::complex<double> >(this); }
