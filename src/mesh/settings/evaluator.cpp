
#include <fstream>

#include "evaluator.h"
#include <numeric>
#include "../../configuration/environment.h"

using namespace espreso;

Evaluator::Evaluator(Property property)
: _name(""), _property(property)
{

}


Evaluator::Evaluator(const std::string &name, Property property)
: _name(name), _property(property)
{

}

Evaluator* Evaluator::create(std::ifstream &is, const Coordinates &coordinates)
{
	Evaluator::Type type;
	is.read(reinterpret_cast<char *>(&type), sizeof(Evaluator::Type));
	Property property;
	is.read(reinterpret_cast<char *>(&property), sizeof(Property));

	switch (type) {
	case Type::DEFAULT: return new Evaluator(property);
	case Type::CONST: return new ConstEvaluator(is, property);
	case Type::COORDINATE: return new CoordinatesEvaluator(is, coordinates, property);
	case Type::TABLE: return new TableEvaluator(is, property);
	case Type::TABLE_INTERPOLATION: return new TableInterpolationEvaluator(is, property);
	case Type::ARRAY: ESINFO(GLOBAL_ERROR) << "Implement loading of Array evaluator"; return NULL;
	default: ESINFO(GLOBAL_ERROR) << "Unknown evaluator type"; return NULL;
	}
}



void Evaluator::store(std::ofstream& os)
{
	Type type = Type::DEFAULT;
	os.write(reinterpret_cast<const char *>(&type), sizeof(Evaluator::Type));
	os.write(reinterpret_cast<const char *>(&_property), sizeof(Property));
}



ConstEvaluator::ConstEvaluator(double value, Property property)
: Evaluator(property), _value(value)
{

}

ConstEvaluator::ConstEvaluator(std::ifstream &is, Property property)
: Evaluator(property)
{
	is.read(reinterpret_cast<char *>(&_value), sizeof(double));
}

void ConstEvaluator::store(std::ofstream& os)
{
	Type type = Type::CONST;
	os.write(reinterpret_cast<const char *>(&type), sizeof(Evaluator::Type));
	os.write(reinterpret_cast<const char *>(&_property), sizeof(Property));
	os.write(reinterpret_cast<const char *>(&_value), sizeof(double));
}



ExpressionEvaluator::ExpressionEvaluator(const std::string &expression, std::vector<std::string> variables, Property property)
: Evaluator(property),
 _expression(environment->OMP_NUM_THREADS, Expression(expression, variables)),
 _values(environment->OMP_NUM_THREADS, std::vector<double>(variables.size()))
{

}

ExpressionEvaluator::ExpressionEvaluator(std::ifstream &is, std::vector<std::string> variables, Property property)
: Evaluator(property)
{
	eslocal size;
	is.read(reinterpret_cast<char *>(&size), sizeof(eslocal));
	char *buffer = new char[size];
	is.read(buffer, size);
	std::string expression(buffer, size);
	delete buffer;
	_expression.resize(environment->OMP_NUM_THREADS, Expression(expression, variables));
	_values.resize(environment->OMP_NUM_THREADS, std::vector<double>(variables.size()));
}

CoordinatesEvaluator::CoordinatesEvaluator(const std::string &expression, const Coordinates &coordinates, Property property)
: ExpressionEvaluator(expression, { "x", "y", "z", "TEMPERATURE", "TIME" }, property), _coordinates(coordinates)
{

}

CoordinatesEvaluator::CoordinatesEvaluator(std::ifstream &is, const Coordinates &coordinates, Property property)
: ExpressionEvaluator(is, { "x", "y", "z", "TEMPERATURE", "TIME" }, property), _coordinates(coordinates)
{

}

void CoordinatesEvaluator::store(std::ofstream& os)
{
	Type type = Type::COORDINATE;
	os.write(reinterpret_cast<const char *>(&type), sizeof(Evaluator::Type));
	os.write(reinterpret_cast<const char *>(&_property), sizeof(Property));
	eslocal size = _expression[0].expression().size();
	os.write(reinterpret_cast<const char *>(&size), sizeof(eslocal));
	os.write(_expression[0].expression().c_str(), _expression[0].expression().size());
}


TableEvaluator::TableEvaluator(
			const std::string &name,
			const std::vector<std::vector<std::vector<double> > > &table,
			const std::vector<TableProperty> &properties,
			const std::vector<std::vector<double> > &axis,
			Property property)
: Evaluator(name, property), _dimension(properties.size()), _table(table), _properties(properties), _axis(axis)
{

}

TableEvaluator::TableEvaluator(std::ifstream &is, Property property)
: Evaluator(property)
{
	is.read(reinterpret_cast<char *>(&_dimension), sizeof(size_t));

	size_t size;
	is.read(reinterpret_cast<char *>(&size), sizeof(size_t));
	_table.resize(size);
	is.read(reinterpret_cast<char *>(&size), sizeof(size_t));
	for (size_t i = 0; i < _table.size(); i++) {
		_table[i].resize(size);
	}
	is.read(reinterpret_cast<char *>(&size), sizeof(size_t));
	for (size_t i = 0; i < _table.size(); i++) {
		for (size_t j = 0; j < _table[i].size(); j++) {
			_table[i][j].resize(size);
			for (size_t k = 0; k < _table[i][j].size(); k++) {
				is.read(reinterpret_cast<char *>(&_table[i][j][k]), sizeof(double));
			}
		}
	}

	is.read(reinterpret_cast<char *>(&size), sizeof(size_t));
	_properties.resize(size);
	for (size_t i = 0; i < _properties.size(); i++) {
		is.read(reinterpret_cast<char *>(&_properties[i]), sizeof(TableProperty));
	}

	is.read(reinterpret_cast<char *>(&size), sizeof(size_t));
	_axis.resize(size);
	is.read(reinterpret_cast<char *>(&size), sizeof(size_t));
	for (size_t i = 0; i < _axis.size(); i++) {
		_axis[i].resize(size);
		for (size_t j = 0; j < _axis[i].size(); j++) {
			is.read(reinterpret_cast<char *>(&_axis[i][j]), sizeof(double));
		}
	}
}

void TableEvaluator::store(std::ofstream& os)
{
	Type type = Type::TABLE;
	os.write(reinterpret_cast<const char *>(&type), sizeof(Evaluator::Type));
	os.write(reinterpret_cast<const char *>(&_property), sizeof(Property));
	size_t size = _dimension;
	os.write(reinterpret_cast<const char *>(&size), sizeof(size_t));

	size = _table.size();
	os.write(reinterpret_cast<const char *>(&size), sizeof(size_t));
	size = _table[0].size();
	os.write(reinterpret_cast<const char *>(&size), sizeof(size_t));
	size = _table[0][0].size();
	os.write(reinterpret_cast<const char *>(&size), sizeof(size_t));
	for (size_t i = 0; i < _table.size(); i++) {
		for (size_t j = 0; j < _table[i].size(); j++) {
			for (size_t k = 0; k < _table[i][j].size(); k++) {
				os.write(reinterpret_cast<const char *>(&_table[i][j][k]), sizeof(double));
			}
		}
	}

	size = _properties.size();
	os.write(reinterpret_cast<const char *>(&size), sizeof(size_t));
	for (size_t i = 0; i < _properties.size(); i++) {
		os.write(reinterpret_cast<const char *>(&_properties[i]), sizeof(TableProperty));
	}

	size = _axis.size();
	os.write(reinterpret_cast<const char *>(&size), sizeof(size_t));
	size = _axis[0].size();
	os.write(reinterpret_cast<const char *>(&size), sizeof(size_t));
	for (size_t i = 0; i < _axis.size(); i++) {
		for (size_t j = 0; j < _axis[i].size(); j++) {
			os.write(reinterpret_cast<const char *>(&_axis[i][j]), sizeof(double));
		}
	}
}

TableInterpolationEvaluator::TableInterpolationEvaluator(
		const std::string &name,
		const std::vector<std::pair<double, double> > &table,
		Property property)
: Evaluator(name, property), _table(table)
{
	if (!_table.size()) {
		ESINFO(GLOBAL_ERROR) << "Interpolation table with zero size.";
	}
}

TableInterpolationEvaluator::TableInterpolationEvaluator(std::ifstream &is, Property property)
: Evaluator(property)
{
	size_t size;
	is.read(reinterpret_cast<char *>(&size), sizeof(size_t));
	_table.resize(size);
	is.read(reinterpret_cast<char *>(_table.data()), 2 * sizeof(size_t) * _table.size());
}

void TableInterpolationEvaluator::store(std::ofstream& os)
{
	Type type = Type::TABLE_INTERPOLATION;
	os.write(reinterpret_cast<const char *>(&type), sizeof(Evaluator::Type));
	os.write(reinterpret_cast<const char *>(&_property), sizeof(Property));

	size_t size = _table.size();
	os.write(reinterpret_cast<const char *>(&size), sizeof(size_t));
	os.write(reinterpret_cast<const char *>(_table.data()), 2 * sizeof(double) * _table.size());
}

ArrayEvaluator::ArrayEvaluator(const std::string &name, std::vector<eslocal> &indices, std::vector<double> &values, eslocal offset, Property property)
: Evaluator(name, property)
{
	std::vector<eslocal> permutation(indices.size());
	std::iota(permutation.begin(), permutation.end(), 0);
	std::sort(permutation.begin(), permutation.end(), [&] (eslocal i, eslocal j) { return indices[i] < indices[j]; });
	_indices.reserve(indices.size());
	_values.reserve(indices.size());
	for (size_t i = 0; i < permutation.size(); i++) {
		_indices.push_back(indices[permutation[i]] - offset);
		_values.push_back(values[permutation[i]]);
	}
}

ArrayEvaluator::ArrayEvaluator(const std::string &name, eslocal size, eslocal *indices, double *values, eslocal offset, Property property)
: Evaluator(name, property)
{
	std::vector<eslocal> permutation(size);
	std::iota(permutation.begin(), permutation.end(), 0);
	std::sort(permutation.begin(), permutation.end(), [&] (eslocal i, eslocal j) { return indices[i] < indices[j]; });
	_indices.reserve(size);
	_values.reserve(size);
	for (size_t i = 0; i < permutation.size(); i++) {
		_indices.push_back(indices[permutation[i]] - offset);
		_values.push_back(values[permutation[i]]);
	}
}

void ArrayEvaluator::addIndex(eslocal index, eslocal value)
{
	size_t offset = std::lower_bound(_indices.begin(), _indices.end(), index) - _indices.begin();
	_indices.insert(_indices.begin() + offset, index);
	_values.insert(_values.begin() + offset, value);
}


void ArrayEvaluator::store(std::ofstream& os)
{
	ESINFO(GLOBAL_ERROR) << "Implement store ArrayEvaluator.";
}





