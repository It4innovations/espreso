
#ifndef SRC_FETI_PRECONDITIONER_WEIGHTFUNCTION_H_
#define SRC_FETI_PRECONDITIONER_WEIGHTFUNCTION_H_

#include "preconditioner.h"

namespace espreso {

template <typename T>
struct WeightFunction: public Preconditioner<T> {
	WeightFunction(FETI<T> &feti);

	void info();
	void update(const step::Step &step);

	void apply(const Vector_Dual<T> &x, Vector_Dual<T> &y);

private:
	using Preconditioner<T>::feti;
	std::vector<Vector_Dense<T> > Btx;
};

}

#endif /* SRC_FETI_PRECONDITIONER_WEIGHTFUNCTION_H_ */
