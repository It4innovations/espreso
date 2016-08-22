
#ifndef SRC_MESH_SETTINGS_EVALUATOR_H_
#define SRC_MESH_SETTINGS_EVALUATOR_H_

#include "esbasis.h"
#include "../structures/coordinates.h"
#include <cilk/cilk.h>
#include <cilk/cilk_api.h>

namespace espreso {

class Coordinates;

class Evaluator {

public:

	virtual Evaluator* copy() const { return new Evaluator(*this); }
	virtual double evaluate(size_t index) const { return 0; }
	virtual ~Evaluator() {};
};

class ConstEvaluator: public Evaluator {

	ConstEvaluator(double value): _value(value) {};
	virtual Evaluator* copy() const { return new ConstEvaluator(*this); }
	double evaluate(size_t index) const { return _value; }

private:
	double _value;
};


/**
 * Evaluator can be called from various CILK workers.
 * Create one instance for each worker in order to avoid race conditions.
 */
class ExpressionEvaluator: public Evaluator {

public:
	ExpressionEvaluator(const std::string &expression, std::vector<std::string> variables, const Coordinates &coordinates)
	: _expression(config::env::CILK_NWORKERS, Expression(expression, variables)), _coordinates(coordinates) {};

	virtual Evaluator* copy() const =0;
	virtual double evaluate(size_t index) const =0;

protected:
	double evaluate(const std::vector<double> &values) const
	{
		return _expression[__cilkrts_get_worker_number()].evaluate(values);
	}

	std::vector<Expression> _expression;
	const Coordinates &_coordinates;
};

class CoordinatesEvaluator: public ExpressionEvaluator {

public:
	CoordinatesEvaluator(const std::string &expression, const Coordinates &coordinates)
	: ExpressionEvaluator(expression, { "x", "y", "z" }, coordinates), _values(config::env::CILK_NWORKERS, std::vector<double>(3)) {};

	virtual Evaluator* copy() const { return new CoordinatesEvaluator(*this); }
	double evaluate(size_t index) const
	{
		_values[__cilkrts_get_worker_number()][0] = _coordinates[index].x;
		_values[__cilkrts_get_worker_number()][1] = _coordinates[index].y;
		_values[__cilkrts_get_worker_number()][2] = _coordinates[index].z;
		return ExpressionEvaluator::evaluate(_values[__cilkrts_get_worker_number()]);
	}

private:
	std::vector<Expression> _expression;
	mutable std::vector<std::vector<double> >_values;
};

}



#endif /* SRC_MESH_SETTINGS_EVALUATOR_H_ */
