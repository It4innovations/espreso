
#ifndef SRC_BASIS_EVALUATOR_PARAMETER_H_
#define SRC_BASIS_EVALUATOR_PARAMETER_H_

#include <string>

namespace espreso {

struct EvaluatorParameter {
	std::string name;
	double value;

	EvaluatorParameter(const std::string &name): name(name), value(0) {}
};

}



#endif /* SRC_BASIS_EVALUATOR_PARAMETER_H_ */
