
#ifndef SRC_PHYSICS_ASSEMBLER_OPERATORS_EXPRESSION_H_
#define SRC_PHYSICS_ASSEMBLER_OPERATORS_EXPRESSION_H_

#include "basis/evaluator/evaluator.h"
#include "analysis/assembler/operator.h"
#include "analysis/assembler/parameter.h"

namespace espreso {

struct ExternalExpression: ActionOperator {
	ExternalExpression(size_t interval, Evaluator *evaluator)
	: evaluator(evaluator)
	{
		isconst = evaluator->parameters.size() == 0;
		action = Action::ASSEMBLE | Action::REASSEMBLE | Action::SOLUTION;
	}

	Evaluator *evaluator;
};

template <size_t nodes, size_t gps, size_t ndim, size_t edim, size_t etype, class Physics, class Setter>
struct ExternalNodeExpression: ExternalExpression, Physics {

	ExternalNodeExpression(size_t interval, Evaluator *evaluator, const Setter &setter)
	: ExternalExpression(interval, evaluator),
	  setter(setter) {}

	Setter setter;

	void simd(typename Physics::Element &element)
	{
		for (size_t n = 0; n < nodes; ++n) {
			for (size_t s = 0; s < SIMD::size; ++s) {
				this->setter(element, n, s, evaluator->evaluate());
			}
		}
	}
};

template <size_t nodes, size_t gps, size_t ndim, size_t edim, size_t etype, class Physics, class Setter>
struct ExternalGPsExpression: ExternalExpression, Physics {

	ExternalGPsExpression(size_t interval, Evaluator *evaluator, const Setter &setter)
	: ExternalExpression(interval, evaluator),
	  setter(setter) {}

	Setter setter;

	void simd(typename Physics::Element &element)
	{
		for (size_t gp = 0; gp < gps; ++gp) {
			for (size_t s = 0; s < SIMD::size; ++s) {
				this->setter(element, gp, s, evaluator->evaluate());
			}
		}
	}
};

}

#endif /* SRC_PHYSICS_ASSEMBLER_OPERATORS_EXPRESSION_H_ */
