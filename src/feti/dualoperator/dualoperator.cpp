
#include "dualoperator.h"
#include "totalfeti.implicit.h"
#include "totalfeti.explicit.h"

#include "esinfo/eslog.hpp"

namespace espreso {

template<typename T>
static DualOperator<T>* _set(FETI<T> *feti)
{
	switch (feti->configuration.method) {
	case FETIConfiguration::METHOD::TOTAL_FETI:
	case FETIConfiguration::METHOD::IMPLICIT_TFETI:
		eslog::info(" = DUAL OPERATOR                                                         IMPLICIT TOTAL FETI = \n");
		return new TotalFETIImplicit<T>(feti);
	case FETIConfiguration::METHOD::EXPLICIT_TFETI:
		eslog::info(" = DUAL OPERATOR                                                         EXPLICIT TOTAL FETI = \n");
		return new TotalFETIExplicit<T>(feti);
	case FETIConfiguration::METHOD::ACCELERATED_TFETI:
		eslog::info(" = DUAL OPERATOR                                                      ACCELERATED TOTAL FETI = \n");
//		return new TotalFETIExplicit<T>(feti);
	case FETIConfiguration::METHOD::HYBRID_FETI:
	default: return nullptr;
	}
}

template <> DualOperator<double>* DualOperator<double>::set(FETI<double> *feti) { return _set<double>(feti); }
template <> DualOperator<std::complex<double> >* DualOperator<std::complex<double> >::set(FETI<std::complex<double> > *feti) { return _set<std::complex<double> >(feti); }


}
