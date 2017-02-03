
#ifndef SRC_ASSEMBLER_STEP_H_
#define SRC_ASSEMBLER_STEP_H_

#include <cstddef>

namespace espreso {

struct Step {
	Step(): load(0) {}

	size_t load;
};

}



#endif /* SRC_ASSEMBLER_STEP_H_ */
