
#ifndef SRC_ANALYSIS_ASSEMBLER_MODULE_PARAMETER_H_
#define SRC_ANALYSIS_ASSEMBLER_MODULE_PARAMETER_H_

#include "basis/evaluator/evaluator.h"

namespace espreso {

template<size_t dimension>
struct ExternalValue {
	std::vector<Evaluator*> evaluator;

	ExternalValue() {};
};

}

#endif /* SRC_ANALYSIS_ASSEMBLER_MODULE_PARAMETER_H_ */
