
#ifndef SRC_INPUT_MESHGENERATOR_PARSERS_EXPRESSION_H_
#define SRC_INPUT_MESHGENERATOR_PARSERS_EXPRESSION_H_

#include <string>
#include "exprtk.hpp"

namespace espreso {

class Expression {

public:

	Expression(const std::string &str);
	Expression(const Expression &other);
	Expression& operator=(const Expression &other);

	double evaluate(double x, double y, double z)
	{
		_x = x;
		_y = y;
		_z = z;
		return _expression.value();
	}

private:
	void parse();

	std::string _str;
	exprtk::symbol_table<double> _symbol_table;
	exprtk::expression<double> _expression;
	double _x, _y, _z;
};

}



#endif /* SRC_INPUT_MESHGENERATOR_PARSERS_EXPRESSION_H_ */
