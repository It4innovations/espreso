
#ifndef SRC_PHYSICS_SOLVER_LOADSTEP_PSEUDOTIMESTEPPINGSOLVER_H_
#define SRC_PHYSICS_SOLVER_LOADSTEP_PSEUDOTIMESTEPPINGSOLVER_H_

#include "loadstepsolver.h"

namespace espreso {

class NonLinearSolverConfiguration;

class PseudoTimeStepping: public LoadStepSolver {

public:
    PseudoTimeStepping(LinearSystem *system, SubStepSolver *subStepSolver, NonLinearSolverConfiguration &configuration, double duration);

    void init(LoadStepSolver *previous);
    void updateStructuralMatrices();

protected:
    void runNextSubstep();

    NonLinearSolverConfiguration &_configuration;
};

}



#endif /* SRC_PHYSICS_SOLVER_LOADSTEP_PSEUDOTIMESTEPPINGSOLVER_H_ */
