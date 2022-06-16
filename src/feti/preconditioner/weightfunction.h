
#ifndef SRC_FETI_PRECONDITIONER_WEIGHTFUNCTION_H_
#define SRC_FETI_PRECONDITIONER_WEIGHTFUNCTION_H_

#include "preconditioner.h"

namespace espreso {

template <typename T>
class WeightFunction: public Preconditioner<T> {
public:
	WeightFunction(FETI<T> *feti);

	void info();
	void update();

	void apply(const Vector_Dual<T> &x, Vector_Dual<T> &y);

	std::vector<Vector_Dense<T> > Btx;
};

}

#endif /* SRC_FETI_PRECONDITIONER_WEIGHTFUNCTION_H_ */