
#ifndef SRC_ASSEMBLER_STEP_H_
#define SRC_ASSEMBLER_STEP_H_

#include <cstddef>

namespace espreso {

struct Step {
	Step(): load(0), time(0), solver(0), subload(0) {}

	size_t load;
	size_t time;
	size_t solver;
	size_t subload;
};

}



#endif /* SRC_ASSEMBLER_STEP_H_ */
