
#include "exprtk.h"
#include "exprtk.hpp"

#include "esinfo/eslog.hpp"

namespace espreso {

struct ExprtkData {
	exprtk::symbol_table<double> symbol_table;
	exprtk::expression<double> expression;

	ExprtkData() {}
	ExprtkData(const ExprtkData &other): symbol_table(other.symbol_table), expression(other.expression) {}
};

template <typename, typename> struct Sequence
{
	Sequence(std::vector<EvaluatorParameter> &variables)
	: variables(variables)
	{

	}

	void push_back(const std::string &parameter)
	{
		variables.push_back({parameter});
	}

	std::vector<EvaluatorParameter> &variables;
};

bool Exprtk::check(const std::string &expression, std::vector<EvaluatorParameter> &variables)
{
	Sequence<std::string, std::string::allocator_type> sequence(variables);

	exprtk::symbol_table<double> symbol_table;
	exprtk::expression<double> expr;
	exprtk::parser<double> parser;

	symbol_table.add_constants();
	exprtk::collect_variables(expression, sequence);
	for (size_t i = 0; i < variables.size(); i++) {
		if (!symbol_table.symbol_exists(variables[i].name)) {
			symbol_table.add_variable(variables[i].name, variables[i].value);
		}
	}

	expr.register_symbol_table(symbol_table);
	return parser.compile(expression, expr);
}

Exprtk::Exprtk(const std::string &expression, std::vector<EvaluatorParameter> &variables)
{
	Sequence<std::string, std::string::allocator_type> sequence(variables);

	_exprtk = new ExprtkData();
	_exprtk->symbol_table.add_constants();

	exprtk::collect_variables(expression, _exprtk->symbol_table, sequence);
	for (size_t i = 0; i < variables.size(); i++) {
		if (!_exprtk->symbol_table.symbol_exists(variables[i].name)) {
			_exprtk->symbol_table.add_variable(variables[i].name, variables[i].value);
		}
	}
	_exprtk->expression.register_symbol_table(_exprtk->symbol_table);

	exprtk::parser<double> parser;
	if (!parser.compile(expression, _exprtk->expression)) {
		eslog::globalerror("invalid expression: '%s'\n", expression.c_str());
	}
}

Exprtk::~Exprtk()
{
	delete _exprtk;
}

double Exprtk::evaluate() const
{
	return _exprtk->expression.value();
}

}


