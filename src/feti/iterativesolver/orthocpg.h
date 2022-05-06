
#ifndef SRC_FETI_ITERATIVESOLVER_ORTHOCPG_H_
#define SRC_FETI_ITERATIVESOLVER_ORTHOCPG_H_

#include "iterativesolver.h"
#include "math/feti/matrix_dual_orthogonal.h"

namespace espreso {

template <typename T>
class OrthogonalizedCPG: public IterativeSolver<T> {
public:
	OrthogonalizedCPG(FETI<T> *feti);

	void info();
	void solve(IterativeSolverInfo &info);

	Vector_Dual<T> l, r, w, x;

	Matrix_Dual_Orthogonal<T> pi, Fpi;
	std::vector<T> wFp, pFp;
};

}


#endif /* SRC_FETI_ITERATIVESOLVER_ORTHOCPG_H_ */
