
#ifndef SRC_ANALYSIS_ASSEMBLER_SUBKERNEL_EXPRESSION_H_
#define SRC_ANALYSIS_ASSEMBLER_SUBKERNEL_EXPRESSION_H_

#include "subkernel.h"

#include <functional>

namespace espreso {

struct ExternalExpression {
	const char* name() const { return expression->value.c_str(); }

	ECFExpression *expression;

	ExternalExpression()
	: expression(nullptr)
	{

	}

	void activate(ECFExpression &expression)
	{
		this->expression = &expression;
	}
};

struct ExternalExpressionVector {
	const char* name() const { return "ExternalExpressionVector"; }

	ECFExpressionVector *expression;

	ExternalExpressionVector()
	: expression(nullptr)
	{

	}

	void activate(ECFExpressionVector &expression)
	{
		this->expression = &expression;
	}
};

struct ExternalEvaluator: SubKernel {
	const char* name() const { return evaluator->expression(); }

	ExternalEvaluator(Evaluator *evaluator)
	: evaluator(evaluator)
	{
		isconst = evaluator->isConst();
		action = Assembler::ASSEMBLE | Assembler::REASSEMBLE | Assembler::ITERATION | Assembler::SOLUTION;
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
struct ExternalGPsExpression: ExternalEvaluator, Physics {

	std::function<void(typename Physics::Element&, size_t&, size_t&, double)> setter;
	int t = omp_get_thread_num();
	double &temperature;
	double &coordinateX, &coordinateY, &coordinateZ;

	ExternalGPsExpression(Evaluator *evaluator, const std::function<void(typename Physics::Element&, size_t&, size_t&, double)> &setter)
	: ExternalEvaluator(evaluator), setter(setter),
	  t(omp_get_thread_num()),
	  temperature(evaluator->getTemperature(t)),
	  coordinateX(evaluator->getCoordinateX(t)), coordinateY(evaluator->getCoordinateY(t)), coordinateZ(evaluator->getCoordinateZ(t))
	{

	}

	void simd(typename Physics::Element &element)
	{
		for (size_t gp = 0; gp < gps; ++gp) {
			for (size_t s = 0; s < SIMD::size; ++s) {
				element.setTemperatureGP(gp, s, temperature);
				element.setCoordinatesGP(gp, s, coordinateX, coordinateY, coordinateZ);
				setter(element, gp, s, this->evaluator->evaluate());
			}
		}
	}
};

template <size_t nodes, class Physics>
struct ExternalNodeExpression: ExternalEvaluator, Physics {

	std::function<void(typename Physics::Element&, size_t&, size_t&, double)> setter;
	int t = omp_get_thread_num();
	double &temperature;
	double &coordinateX, &coordinateY, &coordinateZ;

	ExternalNodeExpression(Evaluator *evaluator, const std::function<void(typename Physics::Element&, size_t&, size_t&, double)> &setter)
	: ExternalEvaluator(evaluator), setter(setter),
	  t(omp_get_thread_num()),
	  temperature(evaluator->getTemperature(t)),
	  coordinateX(evaluator->getCoordinateX(t)), coordinateY(evaluator->getCoordinateY(t)), coordinateZ(evaluator->getCoordinateZ(t))
	{

	}

	void simd(typename Physics::Element &element)
	{
		for (size_t n = 0; n < nodes; ++n) {
			for (size_t s = 0; s < SIMD::size; ++s) {
				element.setTemperatureNode(n, s, temperature);
				element.setCoordinatesNode(n, s, coordinateX, coordinateY, coordinateZ);
				setter(element, n, s, this->evaluator->evaluate());
			}
		}
	}
};

}

#endif /* SRC_ANALYSIS_ASSEMBLER_SUBKERNEL_EXPRESSION_H_ */
