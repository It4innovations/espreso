
#ifndef SRC_PHYSICS_SOLVER_SUBSTEP_SUBSTEPSOLVER_H_
#define SRC_PHYSICS_SOLVER_SUBSTEP_SUBSTEPSOLVER_H_

#include <string>

namespace espreso {

struct LoadStepSolverConfiguration;
class LoadStepSolver;
class LinearSystem;

class SubStepSolver {

	friend class LoadStepSolver;

public:
	SubStepSolver(LinearSystem *system): system(system) {}
	virtual ~SubStepSolver() {}

	virtual void init(SubStepSolver *previous) =0;
	virtual bool hasSameMode(const LoadStepSolverConfiguration &configuration) const =0;
	virtual bool solve(LoadStepSolver &loadStepSolver) =0;

	LinearSystem *system;
};

}



#endif /* SRC_PHYSICS_SOLVER_SUBSTEP_SUBSTEPSOLVER_H_ */
