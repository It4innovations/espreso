
#ifndef SRC_ANALYSIS_LINEARSYSTEM_REGULARIZATION_REGULARIZATION_HEATTRANSFER_H_
#define SRC_ANALYSIS_LINEARSYSTEM_REGULARIZATION_REGULARIZATION_HEATTRANSFER_H_

#include "regularization.h"

namespace espreso {

template <typename T>
struct RegularizationHeatTransfer: Regularization<T> {

	RegularizationHeatTransfer(FETI<T> &feti): Regularization<T>(feti) {}

	void setAnalytic();
	void updateAnalytic();

protected:
	using Regularization<T>::feti;
};

}

#endif /* SRC_ANALYSIS_LINEARSYSTEM_REGULARIZATION_REGULARIZATION_HEATTRANSFER_H_ */
