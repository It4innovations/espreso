
#ifndef SRC_PHYSICS_SOLVER_LOADSTEP_TRANSIENTFIRSTORDERIMPLICITSOLVER_H_
#define SRC_PHYSICS_SOLVER_LOADSTEP_TRANSIENTFIRSTORDERIMPLICITSOLVER_H_

#include "loadstepsolver.h"

namespace espreso {

class TransientFirstOrderImplicitSolverConfiguration;
class Vectors;

class TransientFirstOrderImplicit: public LoadStepSolver {

public:
	TransientFirstOrderImplicit(LinearSystem *system, SubStepSolver *subStepSolver, TransientFirstOrderImplicitSolverConfiguration &configuration, double duration);
	~TransientFirstOrderImplicit();

	void init(LoadStepSolver *previous);
	void updateStructuralMatrices();

protected:
	bool runNextSubstep() override;

	TransientFirstOrderImplicitSolverConfiguration &_configuration;
	double _alpha;
	double _nTimeShift;

	Vectors *U, *dU, *V, *X, *Y, *dTK, *dTM;
};

}



#endif /* SRC_PHYSICS_SOLVER_LOADSTEP_TRANSIENTFIRSTORDERIMPLICITSOLVER_H_ */
