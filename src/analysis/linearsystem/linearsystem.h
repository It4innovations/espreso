
#ifndef SRC_ANALYSIS_LINEARSOLVER_LINEARSOLVER_H_
#define SRC_ANALYSIS_LINEARSOLVER_LINEARSOLVER_H_

#include "analysis/math/matrix_base.h"
#include "analysis/math/vector_base.h"

namespace espreso {

namespace step { struct Step; }

template <typename T>
struct LinearSystemSolver {

	virtual ~LinearSystemSolver() {}

	virtual void set(step::Step &step) =0;
	virtual void update(step::Step &step) =0;
	virtual bool solve(step::Step &step) =0;

	Matrix_Base<T> *A;
	Vector_Base<T> *x, *b, *dirichlet;
};

}

#endif /* SRC_ANALYSIS_LINEARSOLVER_LINEARSOLVER_H_ */
