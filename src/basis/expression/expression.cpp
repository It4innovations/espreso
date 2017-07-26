
#include "expression.h"
#include "../logging/logging.h"

using namespace espreso;

bool Expression::isValid(const std::string &str, std::vector<std::string> variables)
{
	std::vector<double> values(variables.size());

	exprtk::symbol_table<double> symbol_table;
	exprtk::expression<double> expression;
	exprtk::parser<double> parser;

	for (size_t i = 0; i < variables.size(); i++) {
		symbol_table.add_variable(variables[i], values[i]);
	}
	symbol_table.add_constants();
	expression.register_symbol_table(symbol_table);

	return parser.compile(str, expression);
}

void Expression::parse()
{
	if (!isValid(_str, _variables)) {
		ESINFO(GLOBAL_ERROR) << "Invalid expression: '" << _str << "'";
	}

	_values.resize(_variables.size());
	for (size_t i = 0; i < _variables.size(); i++) {
		_symbol_table.add_variable(_variables[i], _values[i]);
	}
	_symbol_table.add_constants();
	_expression.register_symbol_table(_symbol_table);

	exprtk::parser<double> parser;
	parser.compile(_str, _expression);
}

Expression::Expression(const std::string &str, std::vector<std::string> variables): _str(str), _variables(variables)
{
	parse();
}

Expression::Expression(const Expression &other): _str(other._str), _variables(other._variables)
{
	parse();
}

Expression& Expression::operator=(const Expression &other)
{
	if (this != &other) {
		_str = other._str;
		_variables = other._variables;
		parse();
	}
	return *this;
}




