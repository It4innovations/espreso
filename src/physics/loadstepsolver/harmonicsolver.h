
#ifndef SRC_PHYSICS_SOLVER_LOADSTEP_HARMONICSOLVER_H_
#define SRC_PHYSICS_SOLVER_LOADSTEP_HARMONICSOLVER_H_

#include "loadstepsolver.h"

namespace espreso {

struct HarmonicSolverConfiguration;

class HarmonicSolver: public LoadStepSolver {

public:
	HarmonicSolver(LinearSystem *system, SubStepSolver *subStepSolver, HarmonicSolverConfiguration &configuration, double duration);

	void init(LoadStepSolver *previous);
	void updateStructuralMatrices();

protected:
	void runNextSubstep();

	void updateDamping();
	void store();
	void ftt();

	HarmonicSolverConfiguration &_configuration;
	double *_fttRequestedFrequencies;
	double *_fttRequestedFrequenciesEnd;
};

}



#endif /* SRC_PHYSICS_SOLVER_LOADSTEP_HARMONICSOLVER_H_ */
