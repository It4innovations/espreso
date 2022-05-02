
#ifndef SRC_AXFETI_ITERATIVESOLVER_PCPG_H_
#define SRC_AXFETI_ITERATIVESOLVER_PCPG_H_

#include "iterativesolver.h"

namespace espreso {

template <typename T>
class PCPG: public IterativeSolver<T> {
public:
	PCPG(AX_FETI<T> *feti);

	void info();
	void solve(IterativeSolverInfo &info);

	Vector_Dual<T> l, r, w, y, z, p;

	Vector_Dual<T> x, Fp;
};

}

#endif /* SRC_AXFETI_ITERATIVESOLVER_PCPG_H_ */
