
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

	double evaluate() const
	{
		int thread = omp_get_thread_num();
		return _expression[thread]->evaluate();
	}

protected:
	std::vector<Exprtk*> _expression;
};

}


#endif /* SRC_BASIS_EVALUATOR_EXPRESSIONEVALUATOR_H_ */
