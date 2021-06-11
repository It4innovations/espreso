
#ifndef SRC_PHYSICS_SOLVER_LOADSTEP_TRANSIENTSECONDORDERIMPLICITSOLVER_H_
#define SRC_PHYSICS_SOLVER_LOADSTEP_TRANSIENTSECONDORDERIMPLICITSOLVER_H_

#include "loadstepsolver.h"

namespace espreso {

class Vectors;
class TransientSecondOrderImplicitSolverConfiguration;

class TransientSecondOrderImplicit: public LoadStepSolver {

public:
	TransientSecondOrderImplicit(LinearSystem *system, SubStepSolver *subStepSolver, TransientSecondOrderImplicitSolverConfiguration &configuration, double duration);
	~TransientSecondOrderImplicit();

	void init(LoadStepSolver *previous);
	void updateStructuralMatrices();

protected:
	bool runNextSubstep() override;

	void updateConstants();
	void updateDamping();

	TransientSecondOrderImplicitSolverConfiguration &_configuration;
	double _alpha, _delta;
	double _massDamping, _stiffnessDamping;
	double _nTimeShift;
	double _newmarkConsts[8];

	Vectors *U, *dU, *V, *W, *X, *Y, *Z, *dTK, *dTM;
};

}



#endif /* SRC_PHYSICS_SOLVER_LOADSTEP_TRANSIENTSECONDORDERIMPLICITSOLVER_H_ */
