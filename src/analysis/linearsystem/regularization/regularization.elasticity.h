
#ifndef SRC_ANALYSIS_LINEARSYSTEM_REGULARIZATION_REGULARIZATION_ELASTICITY_H_
#define SRC_ANALYSIS_LINEARSYSTEM_REGULARIZATION_REGULARIZATION_ELASTICITY_H_

#include "regularization.h"

namespace espreso {

template <typename T>
struct RegularizationElasticity: Regularization<T> {

	RegularizationElasticity(FETI<T> &feti): Regularization<T>(feti) {}

	void setAnalytic();
	void updateAnalytic();

protected:
	void set2D(esint domain);
	void set3D(esint domain);

	using Regularization<T>::feti;
};

}

#endif /* SRC_ANALYSIS_LINEARSYSTEM_REGULARIZATION_REGULARIZATION_ELASTICITY_H_ */
