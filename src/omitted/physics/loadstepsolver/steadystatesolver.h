
#ifndef SRC_PHYSICS_SOLVER_LOADSTEP_STEADYSTATESOLVER_H_
#define SRC_PHYSICS_SOLVER_LOADSTEP_STEADYSTATESOLVER_H_

#include "loadstepsolver.h"

namespace espreso {

class SteadyStateSolver: public LoadStepSolver {

public:
    SteadyStateSolver(LinearSystem *system, SubStepSolver *subStepSolver, double duration);

    void init(LoadStepSolver *previous);
    void updateStructuralMatrices();

protected:
    void runNextSubstep();
};

}



#endif /* SRC_PHYSICS_SOLVER_LOADSTEP_STEADYSTATESOLVER_H_ */
