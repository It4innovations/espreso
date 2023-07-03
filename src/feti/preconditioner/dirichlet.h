
#ifndef SRC_FETI_PRECONDITIONER_DIRICHLET_H_
#define SRC_FETI_PRECONDITIONER_DIRICHLET_H_

#include "preconditioner.h"

namespace espreso {

template <typename T>
struct Dirichlet: public Preconditioner<T> {
	Dirichlet(FETI<T> *feti);
	~Dirichlet();

	void info();
	void update();

	void apply(const Vector_Dual<T> &x, Vector_Dual<T> &y);

protected:
	void _print();

	std::vector<Vector_Dense<T> > Btx, KBtx;
	std::vector<Matrix_Dense<T> > sc;
	std::vector<DirectSolver<T, Matrix_CSR> > Ksolver;
};

}

#endif /* SRC_FETI_PRECONDITIONER_DIRICHLET_H_ */
