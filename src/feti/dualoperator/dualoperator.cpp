
#include "dualoperator.h"
#include "totalfeti.h"
#include "totalfetiexplicit.h"

namespace espreso {

template<typename T>
static DualOperator<T>* _set(FETI<T> *feti)
{
	switch (feti->configuration.method) {
	case FETIConfiguration::METHOD::TOTAL_FETI:
		eslog::info(" = LINEAR SOLVER                                                                  TOTAL FETI = \n");
		// return new TotalFETI<T>(feti);
		return new TotalFETIExplicit<T>(feti);
	case FETIConfiguration::METHOD::HYBRID_FETI:
	default: return nullptr;
	}
}

template <> DualOperator<double>* DualOperator<double>::set(FETI<double> *feti) { return _set<double>(feti); }
template <> DualOperator<std::complex<double> >* DualOperator<std::complex<double> >::set(FETI<std::complex<double> > *feti) { return _set<std::complex<double> >(feti); }


}
