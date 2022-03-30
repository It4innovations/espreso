
#ifndef SRC_AXFETI_ITERATIVESOLVER_ITERATIVESOLVER_H_
#define SRC_AXFETI_ITERATIVESOLVER_ITERATIVESOLVER_H_

#include "axfeti/feti.h"
#include "math2/feti/vector_dual.h"

namespace espreso {

struct IterativeSolverInfo {
	enum ERROR {

	};

	size_t iterations;
};

template <typename T>
class IterativeSolver {
public:
	static IterativeSolver<T>* set(AX_FETI<T> *feti);

	IterativeSolver(AX_FETI<T> *feti): feti(feti) {}
	virtual ~IterativeSolver() {}

	virtual void info() =0;
	virtual void solve(IterativeSolverInfo &info) =0;

	void reconstructSolution(const Vector_Dual<T> &l, const Vector_Dual<T> &r);

	AX_FETI<T> *feti;
};

}

#endif /* SRC_AXFETI_ITERATIVESOLVER_ITERATIVESOLVER_H_ */
