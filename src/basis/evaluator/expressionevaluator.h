
#ifndef SRC_BASIS_EVALUATOR_EXPRESSIONEVALUATOR_H_
#define SRC_BASIS_EVALUATOR_EXPRESSIONEVALUATOR_H_

#include "evaluator.h"
#include "wrappers/exprtk/exprtk.h"

#include "omp.h"
#include <string>
#include <vector>

namespace espreso {

/**
 * Evaluator can be called from various threads.
 * Create one instance for each worker in order to avoid race conditions.
 */
class ExpressionEvaluator: public Evaluator {

public:
	ExpressionEvaluator(const std::string &expression);
	~ExpressionEvaluator();

	double evaluate(int t) const
	{
		return _expression[t]->evaluate();
	}

	const char* expression() const { return _expr; }

protected:
	const char* _expr;
	std::vector<Exprtk*> _expression;
};

}


#endif /* SRC_BASIS_EVALUATOR_EXPRESSIONEVALUATOR_H_ */
