
#ifndef SRC_ANALYSIS_LINEARSYSTEM_REGULARIZATION_REGULARIZATION_HARMONIC_H_
#define SRC_ANALYSIS_LINEARSYSTEM_REGULARIZATION_REGULARIZATION_HARMONIC_H_

#include "regularization.h"

namespace espreso {

template <typename T>
struct RegularizationHarmonic: Regularization<T> {

	RegularizationHarmonic(FETI<T> &feti): Regularization<T>(feti) {}

	void setAnalytic();
	void updateAnalytic();

protected:
	void set2D(esint domain);
	void set3D(esint domain);

	using Regularization<T>::feti;
	std::vector<Matrix_Dense<T> > NtNNtN;
};

}

#endif /* SRC_ANALYSIS_LINEARSYSTEM_REGULARIZATION_REGULARIZATION_HARMONIC_H_ */
