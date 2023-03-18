
#ifndef SRC_BASIS_EVALUATOR_EVALUATOR_H_
#define SRC_BASIS_EVALUATOR_EVALUATOR_H_

#include "parameter.h"

#include <vector>

namespace espreso {

class Evaluator {
public:
	static Evaluator* create(const std::string &expression);

	virtual ~Evaluator() {};
	virtual double evaluate() const { return 0; }

	std::vector<EvaluatorParameter> parameters;
};

}



#endif /* SRC_BASIS_EVALUATOR_EVALUATOR_H_ */
