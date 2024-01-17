
#include "dualoperator.h"
#include "totalfeti.implicit.h"
#include "totalfeti.explicit.h"
#include "totalfeti.explicit.acc.h"

#include "esinfo/eslog.hpp"

namespace espreso {

template<typename T>
DualOperator<T>* DualOperator<T>::set(FETI<T> &feti, const step::Step &step)
{
	DualOperator<T>* dual = nullptr;
	switch (feti.configuration.dual_operator) {
	case FETIConfiguration::DUAL_OPERATOR::IMPLICIT:
		eslog::info(" = DUAL OPERATOR                                                         IMPLICIT TOTAL FETI = \n");
		dual = new TotalFETIImplicit<T>(feti);
		break;
	case FETIConfiguration::DUAL_OPERATOR::EXPLICIT:
		eslog::info(" = DUAL OPERATOR                                                         EXPLICIT TOTAL FETI = \n");
		dual = new TotalFETIExplicit<T>(feti);
		break;
	case FETIConfiguration::DUAL_OPERATOR::EXPLICIT_GPU:
		if (DirectSolver<Matrix_CSR, T>::provideFactors()) {
			eslog::info(" = DUAL OPERATOR                                                  EXPLICIT TOTAL FETI ON GPU = \n");
			dual = new TotalFETIExplicitAcc<T>(feti);
		} else {
			eslog::globalerror("Linked third party solver does not provide factors that are required for requested dual operator.\n");
		}
		break;
	}
	dual->set(step);
	return dual;
}

template class DualOperator<double>;

}
