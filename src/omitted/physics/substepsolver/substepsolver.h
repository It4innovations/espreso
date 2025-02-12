
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
    SubStepSolver(LinearSystem *system): _system(system) {}
    virtual ~SubStepSolver() {}

    virtual void init(SubStepSolver *previous) =0;
    virtual bool hasSameMode(const LoadStepSolverConfiguration &configuration) const =0;
    virtual void solve(LoadStepSolver &loadStepSolver) =0;

protected:
    LinearSystem *_system;
};

}



#endif /* SRC_PHYSICS_SOLVER_SUBSTEP_SUBSTEPSOLVER_H_ */
