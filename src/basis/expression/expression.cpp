
#include "expression.h"

using namespace espreso;

void Expression::parse()
{
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




