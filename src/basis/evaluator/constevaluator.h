
#ifndef SRC_BASIS_EVALUATOR_CONSTEVALUATOR_H_
#define SRC_BASIS_EVALUATOR_CONSTEVALUATOR_H_

#include "evaluator.h"

namespace espreso {

class ConstEvaluator: public Evaluator {

public:
	ConstEvaluator(double value): _value(value) {}

	double evaluate(int t) const { return _value; }

protected:
	const double _value;
};

}


#endif /* SRC_BASIS_EVALUATOR_CONSTEVALUATOR_H_ */
