
#ifndef SRC_WRAPPERS_EXPRTK_EXPRTK_H_
#define SRC_WRAPPERS_EXPRTK_EXPRTK_H_

#include "basis/evaluator/parameter.h"

#include <string>
#include <vector>

namespace espreso {

struct ExprtkData;

class Exprtk {

public:
	static bool check(const std::string &expression)
	{
		std::vector<EvaluatorParameter> variables;
		return check(expression, variables);
	}
	static bool check(const std::string &expression, std::vector<EvaluatorParameter> &variables);

	Exprtk(const std::string &expression, std::vector<EvaluatorParameter> &variables);

	Exprtk(const Exprtk &other) = delete;
	Exprtk& operator=(const Exprtk &other) = delete;

	~Exprtk();

	double evaluate() const;
protected:
	ExprtkData *_exprtk;
};

}



#endif /* SRC_WRAPPERS_EXPRTK_EXPRTK_H_ */
