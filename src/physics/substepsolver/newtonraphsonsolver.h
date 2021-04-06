
#ifndef SRC_PHYSICS_SOLVER_SUBSTEP_NEWTONRAPHSONSOLVER_H_
#define SRC_PHYSICS_SOLVER_SUBSTEP_NEWTONRAPHSONSOLVER_H_

#include "substepsolver.h"

namespace espreso {

class NonLinearSolverConfiguration;
class Matrix;
class Vectors;

class NewtonRaphson: public SubStepSolver {

public:
	NewtonRaphson(LinearSystem *system, NonLinearSolverConfiguration &configuration);
	~NewtonRaphson();

	void init(SubStepSolver *previous);
	bool hasSameMode(const LoadStepSolverConfiguration &configuration) const;

	void solve(LoadStepSolver &loadStepSolver);

protected:
	NonLinearSolverConfiguration &_configuration;

	Matrix *K;
	Vectors *U, *R, *f, *BC;
	Vectors *lsSolution, *lsRHS, *lsResidual;
};

}



#endif /* SRC_PHYSICS_SOLVER_SUBSTEP_NEWTONRAPHSONSOLVER_H_ */
