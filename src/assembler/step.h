
#ifndef SRC_ASSEMBLER_STEP_H_
#define SRC_ASSEMBLER_STEP_H_

#include <cstddef>

namespace espreso {

struct Step {
	Step(): load(0), iteration(0), solver(0) {}

	size_t load;
	size_t iteration;
	size_t solver;
};

}



#endif /* SRC_ASSEMBLER_STEP_H_ */
