
#include "dirichlet.h"
#include "feti/common/applyB.h"

#include "esinfo/ecfinfo.h"
#include "esinfo/eslog.hpp"
#include "basis/utilities/sysutils.h"

#include <algorithm>

namespace espreso {

template <typename T> static void _print(Dirichlet<T> *dirichlet);

template <typename T>
static void _info(Dirichlet<T> *dirichlet)
{
	if (dirichlet->feti->configuration.exhaustive_info) {
		eslog::info(" = DIRICHLET PRECONDITIONER PROPERTIES                                                       = \n");
		eslog::info(" = ----------------------------------------------------------------------------------------- = \n");
	}
}

template <typename T>
static void _set(Dirichlet<T> *dirichlet)
{
	dirichlet->Btx.resize(dirichlet->feti->K->domains.size());
	dirichlet->KBtx.resize(dirichlet->feti->K->domains.size());
	dirichlet->sc.resize(dirichlet->feti->K->domains.size());

	#pragma omp parallel for
	for (size_t d = 0; d < dirichlet->feti->K->domains.size(); ++d) {
		math::initSolver(dirichlet->feti->K->domains[d]);
		dirichlet->Btx[d].resize(dirichlet->feti->K->domains[d].nrows);
		dirichlet->KBtx[d].resize(dirichlet->feti->K->domains[d].nrows);
		dirichlet->sc[d].shape = Matrix_Shape::UPPER;

		const typename FETI<T>::EqualityConstraints::Domain &L = dirichlet->feti->equalityConstraints->domain[d];
		esint sc_size = L.B1.ncols - *std::min_element(L.B1.cols, L.B1.cols + L.B1.nnz);
		dirichlet->sc[d].resize(sc_size, sc_size);
	}

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

// TODO: combine SC with factorization (the same symbolic factorization?)
template <typename T>
static void _update(Dirichlet<T> *dirichlet)
{
	#pragma omp parallel for
	for (size_t d = 0; d < dirichlet->feti->K->domains.size(); ++d) {
		math::computeSC(dirichlet->feti->K->domains[d], dirichlet->sc[d]);
	}
	_print(dirichlet);
}

template <typename T>
static void _apply(Dirichlet<T> *dirichlet, const Vector_Dual<T> &x, Vector_Dual<T> &y)
{
	#pragma omp parallel for
	for (size_t d = 0; d < dirichlet->feti->K->domains.size(); ++d) {
		applyBt(dirichlet->feti, d, x, dirichlet->Btx[d]);
		Vector_Dense<T> _y, _x;
		_y.vals = dirichlet->KBtx[d].vals + dirichlet->KBtx[d].size - dirichlet->sc[d].nrows;
		_x.vals = dirichlet->Btx[d].vals + dirichlet->Btx[d].size - dirichlet->sc[d].nrows;
		math::apply(_y, T{1}, dirichlet->sc[d], T{0}, _x);
	}
	applyB(dirichlet->feti, dirichlet->KBtx, y);
}

template <typename T>
static void _print(Dirichlet<T> *dirichlet)
{
	if (info::ecf->output.print_matrices) {
		eslog::storedata(" STORE: feti/preconditioner/{Dirichlet}\n");
		for (size_t d = 0; d < dirichlet->feti->K->domains.size(); ++d) {
			math::store(dirichlet->sc[d], utils::filename(utils::debugDirectory(*dirichlet->feti->step) + "/feti/precondition", (std::string("Dirichlet") + std::to_string(d)).c_str()).c_str());
		}
	}
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
