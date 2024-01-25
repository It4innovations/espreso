
#ifndef SRC_FETI_PRECONDITIONER_DIRICHLET_H_
#define SRC_FETI_PRECONDITIONER_DIRICHLET_H_

#include "preconditioner.h"
#include "math/wrappers/math.spsolver.h"

namespace espreso {

template <typename T>
struct Dirichlet: public Preconditioner<T> {
	Dirichlet(FETI<T> &feti);
	~Dirichlet();

	void info();
	void update(const step::Step &step);

	void apply(const Vector_Dual<T> &x, Vector_Dual<T> &y);

protected:
	void _print(const step::Step &step);

	using Preconditioner<T>::feti;
	std::vector<Vector_Dense<T> > Btx, KBtx;
	std::vector<Matrix_Dense<T> > sc;
	std::vector<DirectSparseSolver<T> > Ksolver;
};

}

#endif /* SRC_FETI_PRECONDITIONER_DIRICHLET_H_ */
