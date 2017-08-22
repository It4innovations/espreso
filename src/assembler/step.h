
#ifndef SRC_ASSEMBLER_STEP_H_
#define SRC_ASSEMBLER_STEP_H_

#include <cstddef>

namespace espreso {

struct Step {
	Step(): step(0), substep(0), iteration(0), internalForceReduction(1), currentTime(0), timeStep(0),
			timeIntegrationConstantM(0), timeIntegrationConstantK(1), tangentMatrixCorrection(false) {}

	bool isInitial() const { return step == 0 && substep == 0 && iteration == 0; }

	bool operator==(const Step &other) const { return other.step == step && other.substep == substep && other.iteration == iteration; }

	size_t step;
	size_t substep;
	size_t iteration;

	double internalForceReduction;

	double currentTime;
	double timeStep; // difference between current and previous time
	double timeIntegrationConstantM;
	double timeIntegrationConstantK;
	bool tangentMatrixCorrection;
};

}



#endif /* SRC_ASSEMBLER_STEP_H_ */
