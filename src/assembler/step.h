
#ifndef SRC_ASSEMBLER_STEP_H_
#define SRC_ASSEMBLER_STEP_H_

#include <cstddef>

namespace espreso {

struct Step {
	Step(): step(0), substep(0), iteration(0), currentTime(0), timeStep(0) {}

	size_t step;
	size_t substep;
	size_t iteration;

	double currentTime;
	double timeStep; // difference between current and previous time
};

}



#endif /* SRC_ASSEMBLER_STEP_H_ */
