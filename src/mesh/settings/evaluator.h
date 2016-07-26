
#ifndef SRC_MESH_SETTINGS_EVALUATOR_H_
#define SRC_MESH_SETTINGS_EVALUATOR_H_

#include "esbasis.h"
#include <cilk/cilk.h>
#include <cilk/cilk_api.h>

namespace espreso {

class Evaluator {

public:

	virtual Evaluator* copy() const
	{
		return new Evaluator(*this);
	}

	virtual double evaluate(double x, double y, double z) const
	{
		return 0;
	}

	virtual ~Evaluator() {};
};

class ConstEvaluator: public Evaluator {

	ConstEvaluator(double value): _value(value) {};

	virtual Evaluator* copy() const
	{
		return new ConstEvaluator(*this);
	}

	double evaluate(double x, double y, double z) const
	{
		return _value;
	}

private:
	double _value;
};


/**
 * Evaluator can be called from various CILK workers.
 * Create one instance for each worker in order to avoid race conditions.
 */
class ExpressionEvaluator: public Evaluator {

public:
	ExpressionEvaluator(const std::string &expression): _expression(config::env::CILK_NWORKERS, expression) {};

	virtual Evaluator* copy() const
	{
		return new ExpressionEvaluator(*this);
	}

	double evaluate(double x, double y, double z) const
	{
		return _expression[__cilkrts_get_worker_number()].evaluate(x, y, z);
	}

private:
	std::vector<Expression> _expression;
};

}



#endif /* SRC_MESH_SETTINGS_EVALUATOR_H_ */
