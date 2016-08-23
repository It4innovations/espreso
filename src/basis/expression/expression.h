
#ifndef SRC_INPUT_MESHGENERATOR_PARSERS_EXPRESSION_H_
#define SRC_INPUT_MESHGENERATOR_PARSERS_EXPRESSION_H_

#include <string>
#include "exprtk.hpp"

namespace espreso {

class Expression {

public:

	Expression(const std::string &str, std::vector<std::string> variables);
	Expression(const Expression &other);
	Expression& operator=(const Expression &other);
	double operator()(const std::vector<double> &values) const
	{
		return evaluate(values);
	}
	double evaluate(const std::vector<double> &values) const
	{
		_values = values;
		return _expression.value();
	}

protected:
	void parse();

	std::string _str;
	exprtk::symbol_table<double> _symbol_table;
	exprtk::expression<double> _expression;
	std::vector<std::string> _variables;
	mutable std::vector<double> _values;
};

}



#endif /* SRC_INPUT_MESHGENERATOR_PARSERS_EXPRESSION_H_ */
