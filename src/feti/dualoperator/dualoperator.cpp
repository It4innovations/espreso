
#include "dualoperator.h"
#include "totalfeti.implicit.h"
#include "totalfeti.explicit.h"
#include "totalfeti.explicit.acc.h"

#include "esinfo/eslog.hpp"

namespace espreso {

template class DualOperator<double>;

template<typename T>
DualOperator<T>* DualOperator<T>::set(FETI<T> &feti, const step::Step &step)
{
	DualOperator<T>* dual = nullptr;
	switch (feti.configuration.method) {
	case FETIConfiguration::METHOD::TOTAL_FETI:
	case FETIConfiguration::METHOD::IMPLICIT_TFETI:
		eslog::info(" = DUAL OPERATOR                                                         IMPLICIT TOTAL FETI = \n");
		dual = new TotalFETIImplicit<T>(feti);
		break;
	case FETIConfiguration::METHOD::EXPLICIT_TFETI:
		eslog::info(" = DUAL OPERATOR                                                         EXPLICIT TOTAL FETI = \n");
		dual = new TotalFETIExplicit<T>(feti);
		break;
	case FETIConfiguration::METHOD::ACCELERATED_TFETI:
		if (DirectSolver<T, Matrix_CSR>::provideFactors()) {
			eslog::info(" = DUAL OPERATOR                                                      ACCELERATED TOTAL FETI = \n");
			dual = new TotalFETIExplicitAcc<T>(feti);
		} else {
			eslog::globalerror("Linked third party solver does not provide factors that are required for requested dual operator.\n");
		}
		break;
	case FETIConfiguration::METHOD::HYBRID_FETI:
		break;
	default: break;;
	}
	dual->set(step);
	return dual;
}

}
