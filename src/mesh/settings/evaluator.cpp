
#include "evaluator.h"

using namespace espreso;

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
	case Type::ARRAY: ESINFO(GLOBAL_ERROR) << "Implement loading of Array evaluator"; return NULL;
	default: ESINFO(GLOBAL_ERROR) << "Unknown evaluator type"; return NULL;
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

TableEvaluator::TableEvaluator(std::ifstream &is, Property property): Evaluator(property)
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
