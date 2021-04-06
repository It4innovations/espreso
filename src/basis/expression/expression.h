
#ifndef SRC_INPUT_MESHGENERATOR_PARSERS_EXPRESSION_H_
#define SRC_INPUT_MESHGENERATOR_PARSERS_EXPRESSION_H_

#include <string>
#include <vector>

namespace exprtk {
template <typename T> class symbol_table;
template <typename T> class expression;
}

namespace espreso {

class Expression {

public:

	static bool isValid(const std::string &str, std::vector<std::string> variables);

	Expression(const std::string &str, std::vector<std::string> variables);
	Expression(const Expression &other);
	~Expression();
	Expression& operator=(const Expression &other);
	double operator()(const std::vector<double> &values) const { this->values = values; return evaluate(); }
	double evaluate(const std::vector<double> &values) const { this->values = values; return evaluate(); }
	double evaluate() const;

	std::string expression() const { return _str; }

	mutable std::vector<double> values;
protected:
	void parse();

	std::string _str;
	exprtk::symbol_table<double> *_symbol_table;
	exprtk::expression<double> *_expression;
	std::vector<std::string> _variables;
};

}



#endif /* SRC_INPUT_MESHGENERATOR_PARSERS_EXPRESSION_H_ */
