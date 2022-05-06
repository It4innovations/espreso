
#include "dirichlet.h"
#include "feti/common/applyB.h"

#include "esinfo/eslog.hpp"

#include <algorithm>

namespace espreso {

template <typename T>
static void _info(Dirichlet<T> *dual)
{
	if (dual->feti->configuration.exhaustive_info) {
		eslog::info(" = DIRICHLET PRECONDITIONER PROPERTIES                                                       = \n");
		eslog::info(" = ----------------------------------------------------------------------------------------- = \n");
	}
}

template <typename T>
static void _set(Dirichlet<T> *dirichlet)
{
	dirichlet->Btx.resize(dirichlet->feti->K->domains.size());
	dirichlet->KBtx.resize(dirichlet->feti->K->domains.size());

	#pragma omp parallel for
	for (size_t d = 0; d < dirichlet->feti->K->domains.size(); ++d) {
		math::initSolver(dirichlet->feti->K->domains[d]);
		math::symbolicFactorization(dirichlet->feti->K->domains[d]);
		dirichlet->Btx[d].resize(dirichlet->feti->K->domains[d].nrows);
		dirichlet->KBtx[d].resize(dirichlet->feti->K->domains[d].nrows);
	}

	const auto &dir = dirichlet->feti->K->decomposition->fixedDOFs;
	const auto &gluing = dirichlet->feti->K->decomposition->sharedDOFs;
	dirichlet->surface.resize(gluing.size() + dir.size());
	std::set_union(gluing.begin(), gluing.end(), dir.begin(), dir.end(), dirichlet->surface.begin());

	eslog::checkpointln("FETI: SET DIRICHLET PRECONDITIONER");
}

template <typename T>
static void _free(Dirichlet<T> *dirichlet)
{
	#pragma omp parallel for
	for (size_t d = 0; d < dirichlet->feti->K->domains.size(); ++d) {
		math::freeSolver(dirichlet->feti->K->domains[d]);
	}
}


template <typename T>
static void _update(Dirichlet<T> *dirichlet)
{
	#pragma omp parallel for
	for (size_t d = 0; d < dirichlet->feti->K->domains.size(); ++d) {
		math::numericalFactorization(dirichlet->feti->K->domains[d]);
	}
}

template <typename T>
static void _apply(Dirichlet<T> *dirichlet, const Vector_Dual<T> &x, Vector_Dual<T> &y)
{
	#pragma omp parallel for
	for (size_t d = 0; d < dirichlet->feti->K->domains.size(); ++d) {
		applyBt(dirichlet->feti, d, x, dirichlet->Btx[d]);
		math::solve(dirichlet->feti->K->domains[d], dirichlet->Btx[d], dirichlet->KBtx[d]);
	}
	applyB(dirichlet->feti, dirichlet->KBtx, y);
}


template <> Dirichlet<double>::Dirichlet(FETI<double> *feti): Preconditioner(feti) { _set<double>(this); }
template <> Dirichlet<std::complex<double> >::Dirichlet(FETI<std::complex<double> > *feti): Preconditioner(feti) { _set<std::complex<double> >(this); }

template <> Dirichlet<double>::~Dirichlet() { _free<double>(this); }
template <> Dirichlet<std::complex<double> >::~Dirichlet() { _free<std::complex<double> >(this); }

template <> void Dirichlet<double>::info() { _info<double>(this); }
template <> void Dirichlet<std::complex<double> >::info() { _info<std::complex<double> >(this); }

template <> void Dirichlet<double>::update() { _update<double>(this); }
template <> void Dirichlet<std::complex<double> >::update() { _update<std::complex<double> >(this); }

template <> void Dirichlet<double>::apply(const Vector_Dual<double> &x, Vector_Dual<double> &y) { _apply(this, x, y); }
template <> void Dirichlet<std::complex<double> >::apply(const Vector_Dual<std::complex<double> > &x, Vector_Dual<std::complex<double> > &y) { _apply(this, x, y); }

}
