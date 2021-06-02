
#ifndef SRC_PHYSICS_SOLVER_LOADSTEP_LOADSTEPSOLVER_H_
#define SRC_PHYSICS_SOLVER_LOADSTEP_LOADSTEPSOLVER_H_

#include <string>

namespace espreso {

struct LoadStepSolverConfiguration;
class SubStepSolver;
class LinearSystem;
enum Request: int;

class LoadStepSolver {

public:
	LoadStepSolver(LinearSystem *system, SubStepSolver *subStepSolver, double duration);
	virtual ~LoadStepSolver() {}

	virtual void runNextSubstep() =0;

	virtual void init(LoadStepSolver *previous) =0;
	virtual void updateStructuralMatrices() =0;

	LinearSystem *system;
	SubStepSolver *subStepSolver;
};

}


#endif /* SRC_PHYSICS_SOLVER_LOADSTEP_LOADSTEPSOLVER_H_ */
