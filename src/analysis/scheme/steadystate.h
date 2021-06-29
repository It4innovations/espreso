
#ifndef SRC_ANALYSIS_SCHEME_STEADYSTATE_H_
#define SRC_ANALYSIS_SCHEME_STEADYSTATE_H_

#include "analysis/linearsolver/linearsolver.h"

#include "math2/generalization/vector_base.h"
#include "math2/generalization/matrix_base.h"

namespace espreso {

struct AX_SteadyState {

	void init(DirectSolver<double> *solver);
	void init(FETISolver<double> *solver);
	void init(MultigridSolver<double> *solver);

	Matrix_Base<double> *K;
	Vector_Base<double> *f;
};

}

#endif /* SRC_ANALYSIS_SCHEME_STEADYSTATE_H_ */
