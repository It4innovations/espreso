
#ifndef SRC_PHYSICS_ASSEMBLER_OPERATORS_EXPRESSION_H_
#define SRC_PHYSICS_ASSEMBLER_OPERATORS_EXPRESSION_H_

#include <analysis/assembler/subkernel/operator.h>
#include "basis/evaluator/evaluator.h"
#include "analysis/assembler/parameter.h"

namespace espreso {

template <class Setter>
struct ExternalExpression: ActionOperator {
	const char* name() const { return evaluator->expression(); }

	ExternalExpression(int thread, size_t interval, Evaluator *evaluator, const Setter &setter)
	: evaluator(evaluator), setter(setter), thread(thread)
	{
		isconst = evaluator->parameters[thread].size() == 1; // only dummy
		action = Action::ASSEMBLE | Action::REASSEMBLE | Action::SOLUTION;
	}

	Evaluator *evaluator;
	Setter setter;
	const int thread;

	void setTime(double time, int t)
	{
		evaluator->getTime(t) = time;
	}

	void setFrequency(double frequency, int t)
	{
		evaluator->getFrequency(t) = frequency;
	}
};

template <size_t nodes, size_t gps, size_t ndim, size_t edim, size_t etype, class Physics, class Setter>
struct ExternalNodeExpression: ExternalExpression<Setter>, Physics {
	using ExternalExpression<Setter>::ExternalExpression;

	void simd(typename Physics::Element &element)
	{
		for (size_t n = 0; n < nodes; ++n) {
			for (size_t s = 0; s < SIMD::size; ++s) {
				this->setter(element, n, s, this->evaluator->evaluate());
			}
		}
	}
};

template <size_t nodes, size_t gps, size_t ndim, size_t edim, size_t etype, class Physics, class Setter>
struct ExternalGPsExpression: ExternalExpression<Setter>, Physics {
	using ExternalExpression<Setter>::ExternalExpression;

	void simd(typename Physics::Element &element)
	{
		for (size_t gp = 0; gp < gps; ++gp) {
			for (size_t s = 0; s < SIMD::size; ++s) {
				this->setter(element, gp, s, this->evaluator->evaluate());
			}
		}
	}
};

template <size_t nodes, size_t gps, size_t ndim, size_t edim, size_t etype, class Physics, class Setter> struct ExternalNodeExpressionWithCoordinates;
template <size_t nodes, size_t gps, size_t ndim, size_t edim, size_t etype, class Physics, class Setter> struct ExternalGpsExpressionWithCoordinates;
template <size_t nodes, size_t gps, size_t ndim, size_t edim, size_t etype, class Physics, class Setter> struct ExternalNodeExpressionHeatTransfer;
template <size_t nodes, size_t gps, size_t ndim, size_t edim, size_t etype, class Physics, class Setter> struct ExternalGpsExpressionHeatTransfer;


template <size_t nodes, size_t gps, size_t edim, size_t etype, class Physics, class Setter>
struct ExternalNodeExpressionWithCoordinates<nodes, gps, 2, edim, etype, Physics, Setter>: ExternalExpression<Setter>, Physics {
	using ExternalExpression<Setter>::ExternalExpression;

	void simd(typename Physics::Element &element)
	{
		double &x = this->evaluator->getCoordinateX(this->thread);
		double &y = this->evaluator->getCoordinateY(this->thread);
		for (size_t n = 0; n < nodes; ++n) {
			for (size_t s = 0; s < SIMD::size; ++s) {
				setCoordinates(element, n, s, x, y);
				this->setter(element, n, s, this->evaluator->evaluate());
			}
		}
	}
};

template <size_t nodes, size_t gps, size_t edim, size_t etype, class Physics, class Setter>
struct ExternalNodeExpressionWithCoordinates<nodes, gps, 3, edim, etype, Physics, Setter>: ExternalExpression<Setter>, Physics {
	using ExternalExpression<Setter>::ExternalExpression;

	void simd(typename Physics::Element &element)
	{
		double &x = this->evaluator->getCoordinateX(this->thread);
		double &y = this->evaluator->getCoordinateY(this->thread);
		double &z = this->evaluator->getCoordinateZ(this->thread);
		for (size_t n = 0; n < nodes; ++n) {
			for (size_t s = 0; s < SIMD::size; ++s) {
				setCoordinates(element, n, s, x, y, z);
				this->setter(element, n, s, this->evaluator->evaluate());
			}
		}
	}
};

template <size_t nodes, size_t gps, size_t edim, size_t etype, class Physics, class Setter>
struct ExternalGpsExpressionWithCoordinates<nodes, gps, 2, edim, etype, Physics, Setter>: ExternalExpression<Setter>, Physics {
	using ExternalExpression<Setter>::ExternalExpression;

	void simd(typename Physics::Element &element)
	{
		double &x = this->evaluator->getCoordinateX(this->thread);
		double &y = this->evaluator->getCoordinateY(this->thread);
		for (size_t gp = 0; gp < gps; ++gp) {
			for (size_t s = 0; s < SIMD::size; ++s) {
				setGPCoordinates(element, gp, s, x, y);
				this->setter(element, gp, s, this->evaluator->evaluate());
			}
		}
	}
};

template <size_t nodes, size_t gps, size_t edim, size_t etype, class Physics, class Setter>
struct ExternalGpsExpressionWithCoordinates<nodes, gps, 3, edim, etype, Physics, Setter>: ExternalExpression<Setter>, Physics {
	using ExternalExpression<Setter>::ExternalExpression;

	void simd(typename Physics::Element &element)
	{
		double &x = this->evaluator->getCoordinateX(this->thread);
		double &y = this->evaluator->getCoordinateY(this->thread);
		double &z = this->evaluator->getCoordinateZ(this->thread);
		for (size_t gp = 0; gp < gps; ++gp) {
			for (size_t s = 0; s < SIMD::size; ++s) {
				setGPCoordinates(element, gp, s, x, y, z);
				this->setter(element, gp, s, this->evaluator->evaluate());
			}
		}
	}
};

template <size_t nodes, size_t gps, size_t edim, size_t etype, class Physics, class Setter>
struct ExternalNodeExpressionHeatTransfer<nodes, gps, 2, edim, etype, Physics, Setter>: ExternalExpression<Setter>, Physics {
	using ExternalExpression<Setter>::ExternalExpression;

	void simd(typename Physics::Element &element)
	{
		double &temp = this->evaluator->getTemperature(this->thread);
		double &x = this->evaluator->getCoordinateX(this->thread);
		double &y = this->evaluator->getCoordinateY(this->thread);
		for (size_t n = 0; n < nodes; ++n) {
			for (size_t s = 0; s < SIMD::size; ++s) {
				setTemperature(element, n, s, temp);
				setCoordinates(element, n, s, x, y);
				this->setter(element, n, s, this->evaluator->evaluate());
			}
		}
	}
};

template <size_t nodes, size_t gps, size_t edim, size_t etype, class Physics, class Setter>
struct ExternalNodeExpressionHeatTransfer<nodes, gps, 3, edim, etype, Physics, Setter>: ExternalExpression<Setter>, Physics {
	using ExternalExpression<Setter>::ExternalExpression;

	void simd(typename Physics::Element &element)
	{
		double &temp = this->evaluator->getTemperature(this->thread);
		double &x = this->evaluator->getCoordinateX(this->thread);
		double &y = this->evaluator->getCoordinateY(this->thread);
		double &z = this->evaluator->getCoordinateZ(this->thread);
		for (size_t n = 0; n < nodes; ++n) {
			for (size_t s = 0; s < SIMD::size; ++s) {
				setTemperature(element, n, s, temp);
				setCoordinates(element, n, s, x, y, z);
				this->setter(element, n, s, this->evaluator->evaluate());
			}
		}
	}
};

template <size_t nodes, size_t gps, size_t edim, size_t etype, class Physics, class Setter>
struct ExternalGpsExpressionHeatTransfer<nodes, gps, 2, edim, etype, Physics, Setter>: ExternalExpression<Setter>, Physics {
	using ExternalExpression<Setter>::ExternalExpression;

	void simd(typename Physics::Element &element)
	{
		double &temp = this->evaluator->getTemperature(this->thread);
		double &x = this->evaluator->getCoordinateX(this->thread);
		double &y = this->evaluator->getCoordinateY(this->thread);
		for (size_t gp = 0; gp < gps; ++gp) {
			for (size_t s = 0; s < SIMD::size; ++s) {
				setGPTemperature(element, gp, s, temp);
				setGPCoordinates(element, gp, s, x, y);
				this->setter(element, gp, s, this->evaluator->evaluate());
			}
		}
	}
};

template <size_t nodes, size_t gps, size_t edim, size_t etype, class Physics, class Setter>
struct ExternalGpsExpressionHeatTransfer<nodes, gps, 3, edim, etype, Physics, Setter>: ExternalExpression<Setter>, Physics {
	using ExternalExpression<Setter>::ExternalExpression;

	void simd(typename Physics::Element &element)
	{
		double &temp = this->evaluator->getTemperature(this->thread);
		double &x = this->evaluator->getCoordinateX(this->thread);
		double &y = this->evaluator->getCoordinateY(this->thread);
		double &z = this->evaluator->getCoordinateZ(this->thread);
		for (size_t gp = 0; gp < gps; ++gp) {
			for (size_t s = 0; s < SIMD::size; ++s) {
				setGPTemperature(element, gp, s, temp);
				setGPCoordinates(element, gp, s, x, y, z);
				this->setter(element, gp, s, this->evaluator->evaluate());
			}
		}
	}
};

}

#endif /* SRC_PHYSICS_ASSEMBLER_OPERATORS_EXPRESSION_H_ */
