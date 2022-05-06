
#include "lumped.h"
#include "feti/common/applyB.h"

#include "esinfo/eslog.hpp"

namespace espreso {

template <typename T>
static void _info(Lumped<T> *dual)
{
	if (dual->feti->configuration.exhaustive_info) {
		eslog::info(" = LUMPED PRECONDITIONER PROPERTIES                                                          = \n");

		eslog::info(" = ----------------------------------------------------------------------------------------- = \n");
	}
}

template <typename T>
static void _set(Lumped<T> *lumped)
{
	lumped->Btx.resize(lumped->feti->K->domains.size());
	lumped->KBtx.resize(lumped->feti->K->domains.size());

	#pragma omp parallel for
	for (size_t d = 0; d < lumped->feti->K->domains.size(); ++d) {
		math::commit(lumped->feti->K->domains[d]);
		lumped->Btx[d].resize(lumped->feti->K->domains[d].nrows);
		lumped->KBtx[d].resize(lumped->feti->K->domains[d].nrows);
	}

	eslog::checkpointln("FETI: SET LUMPED PRECONDITIONER");
}

template <typename T>
static void _free(Lumped<T> *lumped)
{
	#pragma omp parallel for
	for (size_t d = 0; d < lumped->feti->K->domains.size(); ++d) {
		math::free(lumped->feti->K->domains[d]);
	}
}

template <typename T>
static void _update(Lumped<T> *lumped)
{

}

template <typename T>
static void _apply(Lumped<T> *lumped, const Vector_Dual<T> &x, Vector_Dual<T> &y)
{
	#pragma omp parallel for
	for (size_t d = 0; d < lumped->feti->K->domains.size(); ++d) {
		applyBt(lumped->feti, d, x, lumped->Btx[d]);
		math::apply(lumped->KBtx[d], T{1}, lumped->feti->K->domains[d], T{0}, lumped->Btx[d]);
	}
	applyB(lumped->feti, lumped->KBtx, y);
}


template <> Lumped<double>::Lumped(FETI<double> *feti): Preconditioner(feti) { _set<double>(this); }
template <> Lumped<std::complex<double> >::Lumped(FETI<std::complex<double> > *feti): Preconditioner(feti) { _set<std::complex<double> >(this); }

template <> Lumped<double>::~Lumped() { _free<double>(this); }
template <> Lumped<std::complex<double> >::~Lumped() { _free<std::complex<double> >(this); }

template <> void Lumped<double>::info() { _info<double>(this); }
template <> void Lumped<std::complex<double> >::info() { _info<std::complex<double> >(this); }

template <> void Lumped<double>::update() { _update<double>(this); }
template <> void Lumped<std::complex<double> >::update() { _update<std::complex<double> >(this); }

template <> void Lumped<double>::apply(const Vector_Dual<double> &x, Vector_Dual<double> &y) { _apply(this, x, y); }
template <> void Lumped<std::complex<double> >::apply(const Vector_Dual<std::complex<double> > &x, Vector_Dual<std::complex<double> > &y) { _apply(this, x, y); }

}
