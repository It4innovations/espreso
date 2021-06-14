
#ifndef SRC_PHYSICS_SOLVER_LOADSTEP_HARMONICSOLVER_H_
#define SRC_PHYSICS_SOLVER_LOADSTEP_HARMONICSOLVER_H_

#include "loadstepsolver.h"

namespace espreso {

struct StructuralMechanicsLoadStepConfiguration;

class HarmonicSolver: public LoadStepSolver {

public:
	HarmonicSolver(LinearSystem *system, SubStepSolver *subStepSolver, StructuralMechanicsLoadStepConfiguration &configuration, double duration);

	void init(LoadStepSolver *previous);
	void updateSystem();
	void updateStructuralMatrices();

protected:
	bool runNextSubstep() override;

	void updateDamping();
	void store();
	void ftt();

	StructuralMechanicsLoadStepConfiguration &_configuration;
	double *_fttRequestedFrequencies;
	double *_fttRequestedFrequenciesEnd;
};

}



#endif /* SRC_PHYSICS_SOLVER_LOADSTEP_HARMONICSOLVER_H_ */
