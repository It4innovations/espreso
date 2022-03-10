
#include "dualoperator.h"
#include "totalfeti.h"

namespace espreso {

template<typename T>
static DualOperator<T>* _set(AX_FETI<T> *feti)
{
	switch (feti->configuration.method) {
	case FETIConfiguration::METHOD::TOTAL_FETI:
		eslog::info(" = LINEAR SOLVER                                                                  TOTAL FETI = \n");
		return new TotalFETI<T>(feti);
		break;
	case FETIConfiguration::METHOD::HYBRID_FETI:
	default: return nullptr;
	}
}

template <> DualOperator<double>* DualOperator<double>::set(AX_FETI<double> *feti) { return _set<double>(feti); }
template <> DualOperator<std::complex<double> >* DualOperator<std::complex<double> >::set(AX_FETI<std::complex<double> > *feti) { return _set<std::complex<double> >(feti); }


}
