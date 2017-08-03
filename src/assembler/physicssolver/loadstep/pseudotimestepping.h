
#ifndef SRC_ASSEMBLER_PHYSICSSOLVER_LOADSTEP_PSEUDOTIMESTEPPING_H_
#define SRC_ASSEMBLER_PHYSICSSOLVER_LOADSTEP_PSEUDOTIMESTEPPING_H_

#include "loadstepsolver.h"

namespace espreso {

class NonLinearSolverBase;

class PseudoTimeStepping: public LoadStepSolver {

public:
	PseudoTimeStepping(TimeStepSolver &timeStepSolver, const NonLinearSolverBase &configuration, double duration);

	Matrices updateStructuralMatrices(Step &step, Matrices matrices);
	Matrices reassembleStructuralMatrices(Step &step, Matrices matrices);

protected:
	void runNextTimeStep(Step &step);
	void processTimeStep(Step &step);

	const NonLinearSolverBase &_configuration;
};

}



#endif /* SRC_ASSEMBLER_PHYSICSSOLVER_LOADSTEP_PSEUDOTIMESTEPPING_H_ */
