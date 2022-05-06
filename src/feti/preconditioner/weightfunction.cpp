
#include "weightfunction.h"
#include "feti/common/applyB.h"

#include "esinfo/eslog.hpp"

namespace espreso {

template <typename T>
static void _info(WeightFunction<T> *dual)
{
	if (dual->feti->configuration.exhaustive_info) {
		eslog::info(" = WEIGHT FUNCTION PRECONDITIONING PROPERTIES                                                = \n");
		eslog::info(" = ----------------------------------------------------------------------------------------- = \n");
	}
}

template <typename T>
static void _set(WeightFunction<T> *wf)
{
	wf->Btx.resize(wf->feti->K->domains.size());

	#pragma omp parallel for
	for (size_t d = 0; d < wf->feti->K->domains.size(); ++d) {
		wf->Btx[d].resize(wf->feti->K->domains[d].nrows);
	}

	eslog::checkpointln("FETI: SET WEIGHT FUNCTION PRECONDITIONER");
}

template <typename T>
static void _update(WeightFunction<T> *wf)
{

}

template <typename T>
static void _apply(WeightFunction<T> *wf, const Vector_Dual<T> &x, Vector_Dual<T> &y)
{
	#pragma omp parallel for
	for (size_t d = 0; d < wf->feti->K->domains.size(); ++d) {
		applyBt(wf->feti, d, x, wf->Btx[d]);
	}
	applyB(wf->feti, wf->Btx, y);
}

template <> WeightFunction<double>::WeightFunction(FETI<double> *feti): Preconditioner(feti) { _set<double>(this); }
template <> WeightFunction<std::complex<double> >::WeightFunction(FETI<std::complex<double> > *feti): Preconditioner(feti) { _set<std::complex<double> >(this); }

template <> void WeightFunction<double>::info() { _info<double>(this); }
template <> void WeightFunction<std::complex<double> >::info() { _info<std::complex<double> >(this); }

template <> void WeightFunction<double>::update() { _update<double>(this); }
template <> void WeightFunction<std::complex<double> >::update() { _update<std::complex<double> >(this); }

template <> void WeightFunction<double>::apply(const Vector_Dual<double> &x, Vector_Dual<double> &y) { _apply(this, x, y); }
template <> void WeightFunction<std::complex<double> >::apply(const Vector_Dual<std::complex<double> > &x, Vector_Dual<std::complex<double> > &y) { _apply(this, x, y); }

}
