
#include "expression.h"

using namespace espreso;

void Expression::parse()
{
	_symbol_table.add_variable("x", _x);
	_symbol_table.add_variable("y", _y);
	_symbol_table.add_variable("z", _z);
	_symbol_table.add_constants();

	_expression.register_symbol_table(_symbol_table);

	exprtk::parser<double> parser;
	parser.compile(_str, _expression);
}

Expression::Expression(const std::string &str): _str(str), _x(0), _y(0), _z(0)
{
	parse();
}

Expression::Expression(const Expression &other): _str(other._str), _x(0), _y(0), _z(0)
{
	parse();
}

Expression& Expression::operator=(const Expression &other)
{
	if (this != &other) {
		_str = other._str;
		_x = 0;
		_y = 0;
		_z = 0;
		parse();
	}
	return *this;
}




