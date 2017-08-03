
#ifndef SRC_ASSEMBLER_PHYSICSSOLVER_LOADSTEP_TRANSIENTFIRSTORDERIMPLICIT_H_
#define SRC_ASSEMBLER_PHYSICSSOLVER_LOADSTEP_TRANSIENTFIRSTORDERIMPLICIT_H_

#include "loadstepsolver.h"

#include <vector>

namespace espreso {

class Solution;
class TransientSolver;

class TransientFirstOrderImplicit: public LoadStepSolver {

public:
	TransientFirstOrderImplicit(TimeStepSolver &timeStepSolver, const TransientSolver &configuration, double duration);

	Matrices updateStructuralMatrices(Step &step, Matrices matrices);
	Matrices reassembleStructuralMatrices(Step &step, Matrices matrices);

protected:
	void initLoadStep(Step &step);
	void runNextTimeStep(Step &step);
	void processTimeStep(Step &step);

	const TransientSolver &_configuration;
	double _alpha;
	bool _timeDependent;

	static size_t loadStep;
	static std::vector<Solution*> solutions;

private:
	enum SolutionIndex { U, dU, V, X, Y };
};

}



#endif /* SRC_ASSEMBLER_PHYSICSSOLVER_LOADSTEP_TRANSIENTFIRSTORDERIMPLICIT_H_ */
