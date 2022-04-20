
#ifndef SRC_AXFETI_ITERATIVESOLVER_ITERATIVESOLVER_H_
#define SRC_AXFETI_ITERATIVESOLVER_ITERATIVESOLVER_H_

#include "axfeti/feti.h"
#include "math/feti/vector_dual.h"

namespace espreso {

struct IterativeSolverInfo {
	enum class ERROR: int {
		OK = 0,
		STAGNATION,
		MAX_ITERATIONS_REACHED,
		INVALID_DATA,
		CONVERGENCE_ERROR,
	};

	ERROR error = ERROR::OK;
	size_t iterations = 0;
	bool converged = false;

	struct Norm {
		struct Dual {
			double absolute, relative, arioli, initial, ksi, criteria;
		} dual;
		double primal;
	} norm;

	struct Time {
		double current, total;
	} time;

	struct Stagnation {
		std::vector<double> buffer;
		int p = 0;
	} stagnation;
};

template <typename T>
class IterativeSolver {
public:
	static IterativeSolver<T>* set(AX_FETI<T> *feti);

	IterativeSolver(AX_FETI<T> *feti): feti(feti) {}
	virtual ~IterativeSolver() {}

	virtual void info() =0;
	virtual void solve(IterativeSolverInfo &info) =0;

	void setInfo(IterativeSolverInfo &info, const FETIConfiguration &configuration, const T &ww);
	void updateInfo(IterativeSolverInfo &info, const FETIConfiguration &configuration, const T &ww, const T &psi, const T &ry);
	void reconstructSolution(const Vector_Dual<T> &l, const Vector_Dual<T> &r);

	AX_FETI<T> *feti;
};

}

#endif /* SRC_AXFETI_ITERATIVESOLVER_ITERATIVESOLVER_H_ */
