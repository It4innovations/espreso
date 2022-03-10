
#include "totalfeti.h"

#include "basis/utilities/sysutils.h"
#include "esinfo/eslog.hpp"
#include "esinfo/ecfinfo.h"

namespace espreso {

template <typename T>
static void _info(TotalFETI<T> *dual)
{
	eslog::info(" = DOMAINS                                                                          %8d = \n", dual->feti->sinfo.domains);
}

template <typename T>
static void _set(TotalFETI<T> *dual)
{
	dual->Kplus.combine(*dual->feti->K, dual->feti->regularization->RegMat);
	dual->Kplus.commit();
	dual->Kplus.symbolicFactorization();
}

template <typename T>
static void _update(TotalFETI<T> *dual)
{
	dual->Kplus.sumCombined(T{1}, *dual->feti->K, dual->feti->regularization->RegMat);
	dual->Kplus.numericalFactorization();

	if (info::ecf->output.print_matrices) {
		eslog::storedata(" STORE: feti/{Kplus}\n");
		math::store(dual->Kplus, utils::filename(utils::debugDirectory(*dual->feti->step) + "/feti", "Kplus").c_str());
	}
}

template <> TotalFETI<double>::TotalFETI(AX_FETI<double> *feti): DualOperator(feti) { _set<double>(this); }
template <> TotalFETI<std::complex<double> >::TotalFETI(AX_FETI<std::complex<double> > *feti): DualOperator(feti) { _set<std::complex<double> >(this); }

template <> void TotalFETI<double>::info() { _info<double>(this); }
template <> void TotalFETI<std::complex<double> >::info() { _info<std::complex<double> >(this); }

template <> void TotalFETI<double>::update() { _update<double>(this); }
template <> void TotalFETI<std::complex<double> >::update() { _update<std::complex<double> >(this); }

}
