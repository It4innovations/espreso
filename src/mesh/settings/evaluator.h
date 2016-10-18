
#ifndef SRC_MESH_SETTINGS_EVALUATOR_H_
#define SRC_MESH_SETTINGS_EVALUATOR_H_

#include "esbasis.h"
#include "property.h"
#include "../structures/coordinates.h"
#include <cilk/cilk.h>
#include <cilk/cilk_api.h>

namespace espreso {

class Coordinates;

class Evaluator {

protected:
	enum class Type: int {
		DEFAULT,
		CONST,
		COORDINATE,
		TABLE,
		ARRAY
	};

public:
	Evaluator(Property property = Property::EMPTY): _name(""), _property(property) {};
	Evaluator(const std::string &name, Property property = Property::EMPTY): _name(name), _property(property) {};

	static Evaluator* create(std::ifstream &is, const Coordinates &coordinates);

	virtual Evaluator* copy() const { return new Evaluator(*this); }
	virtual double evaluate(size_t index, size_t timeStep = 1, double temperature = 0, double pressure = 0, double velocity = 0) const { return 0; }
	virtual double evaluate(const Point &p) const { return 0; }
	virtual const std::string& name() const { return _name; }
	virtual ~Evaluator() {};

	Property& property() { return _property; }
	const Property& property() const { return _property; }

	virtual void store(std::ofstream& os)
	{
		Type type = Type::DEFAULT;
		os.write(reinterpret_cast<const char *>(&type), sizeof(Evaluator::Type));
		os.write(reinterpret_cast<const char *>(&_property), sizeof(Property));
	}

protected:
	std::string _name;
	Property _property;
};

class ConstEvaluator: public Evaluator {

public:
	ConstEvaluator(double value, Property property = Property::EMPTY): Evaluator(property), _value(value) {};
	ConstEvaluator(std::ifstream &is, Property property): Evaluator(property) { is.read(reinterpret_cast<char *>(&_value), sizeof(double)); }

	virtual Evaluator* copy() const { return new ConstEvaluator(*this); }
	double evaluate(size_t index, size_t timeStep, double temperature, double pressure, double velocity) const { return _value; }
	double evaluate(const Point &p) const { return _value; }

	virtual void store(std::ofstream& os)
	{
		Type type = Type::CONST;
		os.write(reinterpret_cast<const char *>(&type), sizeof(Evaluator::Type));
		os.write(reinterpret_cast<const char *>(&_property), sizeof(Property));
		os.write(reinterpret_cast<const char *>(&_value), sizeof(double));
	}

private:
	double _value;
};


/**
 * Evaluator can be called from various CILK workers.
 * Create one instance for each worker in order to avoid race conditions.
 */
class ExpressionEvaluator: public Evaluator {

protected:
	ExpressionEvaluator(const std::string &expression, std::vector<std::string> variables, Property property = Property::EMPTY)
	: Evaluator(property),
	 _expression(config::env::CILK_NWORKERS, Expression(expression, variables)),
	 _values(config::env::CILK_NWORKERS, std::vector<double>(variables.size())) {};

	ExpressionEvaluator(std::ifstream &is, std::vector<std::string> variables, Property property)
	: Evaluator(property)
	{
		eslocal size;
		is.read(reinterpret_cast<char *>(&size), sizeof(eslocal));
		char *buffer = new char[size];
		is.read(buffer, size);
		std::string expression(buffer, size);
		delete buffer;
		_expression.resize(config::env::CILK_NWORKERS, Expression(expression, variables));
		_values.resize(config::env::CILK_NWORKERS, std::vector<double>(variables.size()));
	}

	virtual Evaluator* copy() const =0;
	virtual double evaluate(size_t index, size_t timeStep, double temperature, double pressure, double velocity) const =0;
	virtual double evaluate(const Point &p) const =0;

	std::vector<Expression> _expression;
	mutable std::vector<std::vector<double> >_values;
};

class CoordinatesEvaluator: public ExpressionEvaluator {

public:
	CoordinatesEvaluator(const std::string &expression, const Coordinates &coordinates, Property property = Property::EMPTY)
	: ExpressionEvaluator(expression, { "x", "y", "z" }, property), _coordinates(coordinates) {};

	CoordinatesEvaluator(std::ifstream &is, const Coordinates &coordinates, Property property)
	: ExpressionEvaluator(is, { "x", "y", "z" }, property), _coordinates(coordinates) {};

	virtual Evaluator* copy() const { return new CoordinatesEvaluator(*this); }
	double evaluate(size_t index, size_t timeStep, double temperature, double pressure, double velocity) const
	{
		_values[__cilkrts_get_worker_number()][0] = _coordinates[index].x;
		_values[__cilkrts_get_worker_number()][1] = _coordinates[index].y;
		_values[__cilkrts_get_worker_number()][2] = _coordinates[index].z;
		return _expression[__cilkrts_get_worker_number()].evaluate(_values[__cilkrts_get_worker_number()]);
	}

	double evaluate(const Point &p) const
	{
		_values[__cilkrts_get_worker_number()][0] = p.x;
		_values[__cilkrts_get_worker_number()][1] = p.y;
		_values[__cilkrts_get_worker_number()][2] = p.z;
		return _expression[__cilkrts_get_worker_number()].evaluate(_values[__cilkrts_get_worker_number()]);
	}

	virtual void store(std::ofstream& os)
	{
		Type type = Type::COORDINATE;
		os.write(reinterpret_cast<const char *>(&type), sizeof(Evaluator::Type));
		os.write(reinterpret_cast<const char *>(&_property), sizeof(Property));
		eslocal size = _expression[0].expression().size();
		os.write(reinterpret_cast<const char *>(&size), sizeof(eslocal));
		os.write(_expression[0].expression().c_str(), _expression[0].expression().size());
	}

private:
	const Coordinates &_coordinates;
};

class TableEvaluator: public Evaluator {

public:
	enum class TableProperty {
		TIME,
		TEMPERATURE,
		PRESSURE,
		VELOCITY
	};

	TableEvaluator(
			const std::string &name,
			const std::vector<std::vector<std::vector<double> > > &table,
			const std::vector<TableProperty> &properties,
			const std::vector<std::vector<double> > &axis,
			Property property = Property::EMPTY)
	: Evaluator(name, property), _dimension(properties.size()), _table(table), _properties(properties), _axis(axis) {};

	TableEvaluator(std::ifstream &is, Property property);

	virtual Evaluator* copy() const { return new TableEvaluator(*this); }
	virtual double evaluate(size_t index, size_t timeStep, double temperature, double pressure, double velocity) const
	{
		std::vector<size_t> cell(_dimension);

		for (size_t i = 0; i < _dimension; i++) {
			switch (_properties[i]) {
			case TableProperty::TIME:
				cell[i] = std::find(_axis[i].begin(), _axis[i].end(), timeStep) - _axis[i].begin();
				break;
			case TableProperty::TEMPERATURE:
				cell[i] = std::find(_axis[i].begin(), _axis[i].end(), temperature) - _axis[i].begin();
				break;
			case TableProperty::PRESSURE:
				cell[i] = std::find(_axis[i].begin(), _axis[i].end(), pressure) - _axis[i].begin();
				break;
			case TableProperty::VELOCITY:
				cell[i] = std::find(_axis[i].begin(), _axis[i].end(), velocity) - _axis[i].begin();
				break;
			}
		}
		return _table[_dimension > 0 ? cell[0] : 0][_dimension > 1 ? cell[1] : 0][_dimension > 2 ? cell[2] : 0];
	}

	virtual void store(std::ofstream& os);

protected:
	size_t _dimension;
	std::vector<std::vector<std::vector<double> > > _table;
	std::vector<TableProperty> _properties;
	std::vector<std::vector<double> > _axis;
};

namespace input { class API; }

class ArrayEvaluator: public Evaluator {

	friend class input::API;
public:
	ArrayEvaluator(const std::string &name, std::vector<eslocal> &indices, std::vector<double> &values, eslocal offset, Property property = Property::EMPTY)
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

	ArrayEvaluator(const std::string &name, eslocal size, eslocal *indices, double *values, eslocal offset, Property property = Property::EMPTY)
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


	virtual Evaluator* copy() const { return new ArrayEvaluator(*this); }
	virtual double evaluate(size_t index, size_t timeStep, double temperature, double pressure, double velocity) const
	{
		auto it = std::lower_bound(_indices.begin(), _indices.end(), index);
		if (it != _indices.end() && *it == index) {
			return _values[it - _indices.begin()];
		} else {
			ESINFO(ERROR) << "Array evaluator has no specified value for index '" << index << "'";
			return 0;
		}
	}

	virtual void store(std::ofstream& os)
	{
		Type type = Type::ARRAY;
		ESINFO(GLOBAL_ERROR) << "implement store Array evaluator";
	}

protected:
	std::vector<eslocal> _indices;
	std::vector<double> _values;

};

}



#endif /* SRC_MESH_SETTINGS_EVALUATOR_H_ */
