
#ifndef SRC_ANALYSIS_ASSEMBLER_SUBKERNEL_EXPRESSION_H_
#define SRC_ANALYSIS_ASSEMBLER_SUBKERNEL_EXPRESSION_H_

#include "subkernels.h"

#include <functional>

namespace espreso {

template <size_t gps, class Physics>
struct ExternalGPsExpressionKernel: ExternalExpressionKernel, Physics {

	std::function<void(typename Physics::Element&, size_t&, size_t&, double)> setter;

	ExternalGPsExpressionKernel(Evaluator *evaluator, const std::function<void(typename Physics::Element&, size_t&, size_t&, double)> &setter)
	: ExternalExpressionKernel(evaluator), setter(setter)
	{

	}

	void simd(typename Physics::Element &element)
	{
		for (size_t gp = 0; gp < gps; ++gp) {
			for (size_t s = 0; s < SIMD::size; ++s) {
				setter(element, gp, s, this->evaluator->evaluate());
			}
		}
	}
};

}

#endif /* SRC_ANALYSIS_ASSEMBLER_SUBKERNEL_EXPRESSION_H_ */
