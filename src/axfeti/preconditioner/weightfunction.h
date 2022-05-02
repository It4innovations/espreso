
#ifndef SRC_AXFETI_PRECONDITIONER_WEIGHTFUNCTION_H_
#define SRC_AXFETI_PRECONDITIONER_WEIGHTFUNCTION_H_

#include "preconditioner.h"

namespace espreso {

template <typename T>
class WeightFunction: public Preconditioner<T> {
public:
	WeightFunction(AX_FETI<T> *feti): Preconditioner<T>(feti) {}

	void info() {}
	void update() {}
};

}

#endif /* SRC_AXFETI_PRECONDITIONER_WEIGHTFUNCTION_H_ */
