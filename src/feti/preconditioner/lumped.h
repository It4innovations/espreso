
#ifndef SRC_FETI_PRECONDITIONER_LUMPED_H_
#define SRC_FETI_PRECONDITIONER_LUMPED_H_

#include "preconditioner.h"
#include "math/wrappers/math.spblas.h"

namespace espreso {

template <typename T>
struct Lumped: public Preconditioner<T> {
	Lumped(FETI<T> &feti);
	~Lumped();

	void info();
	void update(const step::Step &step);

	void apply(const Vector_Dual<T> &x, Vector_Dual<T> &y);

protected:
	using Preconditioner<T>::feti;
	std::vector<Vector_Dense<T> > Btx, KBtx;
	std::vector<SpBLAS<Matrix_CSR, T> > KSpBlas;
};

}

#endif /* SRC_FETI_PRECONDITIONER_LUMPED_H_ */
