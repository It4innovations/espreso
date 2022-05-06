
#ifndef SRC_FETI_PRECONDITIONER_WEIGHTFUNCTION_H_
#define SRC_FETI_PRECONDITIONER_WEIGHTFUNCTION_H_

#include "preconditioner.h"

namespace espreso {

template <typename T>
class WeightFunction: public Preconditioner<T> {
public:
	WeightFunction(FETI<T> *feti): Preconditioner<T>(feti) {}

	void info() {}
	void update() {}
};

}

#endif /* SRC_FETI_PRECONDITIONER_WEIGHTFUNCTION_H_ */
