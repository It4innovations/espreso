
#ifndef SRC_PHYSICS_ASSEMBLER_OPERATORS_EXPRESSION_H_
#define SRC_PHYSICS_ASSEMBLER_OPERATORS_EXPRESSION_H_

#include "basis/evaluator/evaluator.h"
#include "basis/expression/variable.h"
#include "analysis/assembler/operator.h"
#include "analysis/assembler/parameter.h"

namespace espreso {

template <class Setter>
struct ExternalExpression: ActionOperator {
	ExternalExpression(size_t interval, Evaluator *evaluator, const Setter &setter)
	: evaluator(evaluator),
	  params(evaluator->params),
	  setter(setter)
	{
		for (size_t i = 0; i < params.general.size(); ++i) {
			params.general[i].variable->set(interval, params.general[i]);
		}
		isconst = evaluator->variables.size() == 0;
		action = Action::ASSEMBLE | Action::REASSEMBLE | Action::SOLUTION;
	}

	Evaluator *evaluator;
	Evaluator::Params params;
	Setter setter;

	void move(int n)
	{
		for (size_t i = 0; i < params.general.size(); ++i) {
			params.general[i].val += n * params.general[i].increment;
		}
	}
};

template <size_t nodes, size_t gps, size_t ndim, size_t edim, size_t etype, class Physics, class Setter>
struct ExternalNodeExpression: ExternalExpression<Setter>, Physics {
	using ExternalExpression<Setter>::ExternalExpression;

	void simd(typename Physics::Element &element)
	{
		double results[SIMD::size * nodes];
		this->evaluator->evalVector(SIMD::size * nodes, this->params, results);
		for (size_t n = 0; n < nodes; ++n) {
			for (size_t s = 0; s < SIMD::size; ++s) {
				this->setter(element, n, s, results[nodes * s + n]); // TODO: check correct order
			}
		}
		this->move(SIMD::size * nodes);
	}
};

template <size_t nodes, size_t gps, size_t ndim, size_t edim, size_t etype, class Physics, class Setter>
struct ExternalGPsExpression: ExternalExpression<Setter>, Physics {
	using ExternalExpression<Setter>::ExternalExpression;

	void simd(typename Physics::Element &element)
	{
		double results[SIMD::size * gps];
		this->evaluator->evalVector(SIMD::size * gps, this->params, results);
		for (size_t gp = 0; gp < gps; ++gp) {
			for (size_t s = 0; s < SIMD::size; ++s) {
				this->setter(element, gp, s, results[gps * s + gp]); // TODO: check correct order
			}
		}
		this->move(SIMD::size * gps);
	}
};

}

#endif /* SRC_PHYSICS_ASSEMBLER_OPERATORS_EXPRESSION_H_ */
