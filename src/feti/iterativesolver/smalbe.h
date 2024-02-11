
#ifndef SRC_FETI_ITERATIVESOLVER_SMALSE_H_
#define SRC_FETI_ITERATIVESOLVER_SMALSE_H_

#include "iterativesolver.h"

namespace espreso {

// Semi-Monotonic Augmented Lagrangian algorithm for Bound and Equality Constraints

template <typename T>
class SMALBE: public IterativeSolver<T> {
public:
	SMALBE(FETI<T> &feti);

	void info();
	void solve(const step::Step &step, IterativeSolverInfo &info);

	using IterativeSolver<T>::feti;
	Vector_Dual<T> l, r, w, p;
	Vector_Dual<T> x, Fp;
};

}


#endif /* SRC_FETI_ITERATIVESOLVER_SMALSE_H_ */