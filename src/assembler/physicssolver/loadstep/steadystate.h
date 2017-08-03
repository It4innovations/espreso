
#ifndef SRC_ASSEMBLER_PHYSICSSOLVER_LOADSTEP_STEADYSTATE_H_
#define SRC_ASSEMBLER_PHYSICSSOLVER_LOADSTEP_STEADYSTATE_H_

#include "loadstepsolver.h"

namespace espreso {

class SteadyStateSolver: public LoadStepSolver {

public:
	SteadyStateSolver(TimeStepSolver &timeStepSolver, double duration);

	Matrices updateStructuralMatrices(Step &step, Matrices matrices);
	Matrices reassembleStructuralMatrices(Step &step, Matrices matrices);

protected:
	void runNextTimeStep(Step &step);
	void processTimeStep(Step &step);
};

}



#endif /* SRC_ASSEMBLER_PHYSICSSOLVER_LOADSTEP_STEADYSTATE_H_ */
