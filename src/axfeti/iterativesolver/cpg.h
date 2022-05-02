
#ifndef SRC_AXFETI_ITERATIVESOLVER_CPG_H_
#define SRC_AXFETI_ITERATIVESOLVER_CPG_H_

#include "iterativesolver.h"

namespace espreso {

template <typename T>
class CPG: public IterativeSolver<T> {
public:
	CPG(AX_FETI<T> *feti);

	void info();
	void solve(IterativeSolverInfo &info);

	Vector_Dual<T> l, r, w, p;

	Vector_Dual<T> x, Fp;
};

}


#endif /* SRC_AXFETI_ITERATIVESOLVER_CPG_H_ */
