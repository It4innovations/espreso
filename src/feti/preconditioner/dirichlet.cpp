
#include "dirichlet.h"

#include "esinfo/eslog.hpp"

namespace espreso {

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

template <typename T>
static void _apply(Dirichlet<T> *dual)
{

}


template <> Dirichlet<double>::Dirichlet(FETI<double> *feti): Preconditioner(feti) { _set<double>(this); }
template <> Dirichlet<std::complex<double> >::Dirichlet(FETI<std::complex<double> > *feti): Preconditioner(feti) { _set<std::complex<double> >(this); }

template <> void Dirichlet<double>::info() { _info<double>(this); }
template <> void Dirichlet<std::complex<double> >::info() { _info<std::complex<double> >(this); }

template <> void Dirichlet<double>::update() { _update<double>(this); }
template <> void Dirichlet<std::complex<double> >::update() { _update<std::complex<double> >(this); }

template <> void Dirichlet<double>::apply(const Vector_Dual<double> &x, Vector_Dual<double> &y) { _apply(this); }
template <> void Dirichlet<std::complex<double> >::apply(const Vector_Dual<std::complex<double> > &x, Vector_Dual<std::complex<double> > &y) { _apply(this); }

}
