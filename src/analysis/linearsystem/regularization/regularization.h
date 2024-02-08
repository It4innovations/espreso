
#ifndef SRC_ANALYSIS_LINEARSYSTEM_REGULARIZATION_REGULARIZATION_H_
#define SRC_ANALYSIS_LINEARSYSTEM_REGULARIZATION_REGULARIZATION_H_

#include "feti/feti.h"

namespace espreso {

template <typename T>
struct Regularization {

	static void set(const step::Step &step, FETI<T> &feti);
	static void update(const step::Step &step, FETI<T> &feti);
};

}

#endif /* SRC_ANALYSIS_LINEARSYSTEM_REGULARIZATION_REGULARIZATION_H_ */
