
#ifndef SRC_ANALYSIS_MODE_LINEAR_H_
#define SRC_ANALYSIS_MODE_LINEAR_H_

#include "analysis/linearsolver/linearsolver.h"

namespace espreso {

struct AX_Linear {

	void init(DirectSolver<double> *solver);
	void init(FETISolver<double> *solver);
	void init(MultigridSolver<double> *solver);

private:
	LinearSolver *solver;
};

}

#endif /* SRC_ANALYSIS_MODE_LINEAR_H_ */
