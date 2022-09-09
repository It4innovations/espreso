
#include "expression.h"
#include "esinfo/eslog.hpp"

#include "exprtk.hpp"

#include <sstream>

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

bool Expression::collectVariables(const std::string &str, std::vector<std::string> &variables)
{
	return exprtk::collect_variables(str, variables);
}

void Expression::parse()
{
	if (!isValid(_str, _variables)) {
		eslog::globalerror("Invalid expression: '%s'\n", _str.c_str());
	}

	values.resize(_variables.size());
	for (size_t i = 0; i < _variables.size(); i++) {
		_symbol_table->add_variable(_variables[i], values[i]);
	}
	_symbol_table->add_constants();
	_expression->register_symbol_table(*_symbol_table);

	exprtk::parser<double> parser;
	parser.compile(_str, *_expression);
}

Expression::Expression(const std::string &str, std::vector<std::string> variables): _str(str), _variables(variables)
{
	_symbol_table = new exprtk::symbol_table<double>();
	_expression = new exprtk::expression<double>();
	parse();
}

Expression::Expression(const Expression &other): _str(other._str), _variables(other._variables)
{
	_symbol_table = new exprtk::symbol_table<double>();
	_expression = new exprtk::expression<double>();
	parse();
}

Expression::~Expression()
{
	delete _symbol_table;
	delete _expression;
}

Expression& Expression::operator=(const Expression &other)
{
	if (this != &other) {
		_str = other._str;
		_variables = other._variables;
		delete _symbol_table;
		delete _expression;
		_symbol_table = new exprtk::symbol_table<double>();
		_expression = new exprtk::expression<double>();
		parse();
	}
	return *this;
}

double Expression::evaluate() const
{
	return _expression->value();
}



