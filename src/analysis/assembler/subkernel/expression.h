
#ifndef SRC_ANALYSIS_ASSEMBLER_SUBKERNEL_EXPRESSION_H_
#define SRC_ANALYSIS_ASSEMBLER_SUBKERNEL_EXPRESSION_H_

#include "subkernel.h"

#include <functional>

namespace espreso {

struct ExternalExpression: SubKernel {
	const char* name() const { return evaluator->expression(); }

	ExternalExpression(Evaluator *evaluator)
	: evaluator(evaluator)
	{
		isconst = evaluator->isConst();
		action = Assembler::ASSEMBLE | Assembler::REASSEMBLE | Assembler::SOLUTION;
	}

	Evaluator *evaluator;

	void setTime(double time, int t)
	{
		evaluator->getTime(t) = time;
	}

	void setFrequency(double frequency, int t)
	{
		evaluator->getFrequency(t) = frequency;
	}
};

template <size_t gps, class Physics>
struct ExternalGPsExpression: ExternalExpression, Physics {

	std::function<void(typename Physics::Element&, size_t&, size_t&, double)> setter;

	ExternalGPsExpression(Evaluator *evaluator, const std::function<void(typename Physics::Element&, size_t&, size_t&, double)> &setter)
	: ExternalExpression(evaluator), setter(setter)
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
