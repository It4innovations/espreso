
#ifndef SRC_PHYSICS_SOLVER_SUBSTEP_LINEARSUBSTEPSOLVER_H_
#define SRC_PHYSICS_SOLVER_SUBSTEP_LINEARSUBSTEPSOLVER_H_

#include "substepsolver.h"

namespace espreso {

class LinearSubStep: public SubStepSolver {

public:
	LinearSubStep(LinearSystem *system);

	void init(SubStepSolver *previous);
	bool hasSameMode(const LoadStepSolverConfiguration &configuration) const;

	void solve(LoadStepSolver &loadStepSolver);
};

}

#endif /* SRC_PHYSICS_SOLVER_SUBSTEP_LINEARSUBSTEPSOLVER_H_ */
