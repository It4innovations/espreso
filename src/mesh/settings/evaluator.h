
#ifndef SRC_MESH_SETTINGS_EVALUATOR_H_
#define SRC_MESH_SETTINGS_EVALUATOR_H_

#include "property.h"
#include "../structures/coordinates.h"

namespace espreso {

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
	Evaluator(Property property = Property::EMPTY);
	Evaluator(const std::string &name, Property property = Property::EMPTY);

	static Evaluator* create(std::ifstream &is, const Coordinates &coordinates);

	virtual Evaluator* copy() const { return new Evaluator(*this); }
	virtual double evaluate(eslocal index, size_t timeStep = 1, double temperature = 0, double pressure = 0, double velocity = 0) const { return 0; }
	virtual double evaluate(const Point &p) const { return 0; }
	virtual const std::string& name() const { return _name; }
	virtual ~Evaluator() {};

	Property& property() { return _property; }
	const Property& property() const { return _property; }

	virtual void store(std::ofstream& os);

protected:
	std::string _name;
	Property _property;
};

class ConstEvaluator: public Evaluator {

public:
	ConstEvaluator(double value, Property property = Property::EMPTY);
	ConstEvaluator(std::ifstream &is, Property property);

	virtual Evaluator* copy() const { return new ConstEvaluator(*this); }
	double inline evaluate(eslocal index, size_t timeStep, double temperature, double pressure, double velocity) const { return _value; }
	double inline evaluate(const Point &p) const { return _value; }

	virtual void store(std::ofstream& os);

private:
	double _value;
};


/**
 * Evaluator can be called from various threads.
 * Create one instance for each worker in order to avoid race conditions.
 */
class ExpressionEvaluator: public Evaluator {

protected:
	ExpressionEvaluator(const std::string &expression, std::vector<std::string> variables, Property property = Property::EMPTY);
	ExpressionEvaluator(std::ifstream &is, std::vector<std::string> variables, Property property);

	virtual Evaluator* copy() const =0;
	virtual inline double evaluate(eslocal index, size_t timeStep, double temperature, double pressure, double velocity) const =0;
	virtual inline double evaluate(const Point &p) const =0;

	std::vector<Expression> _expression;
	mutable std::vector<std::vector<double> >_values;
};

class CoordinatesEvaluator: public ExpressionEvaluator {

public:
	CoordinatesEvaluator(const std::string &expression, const Coordinates &coordinates, Property property = Property::EMPTY);
	CoordinatesEvaluator(std::ifstream &is, const Coordinates &coordinates, Property property);

	virtual Evaluator* copy() const { return new CoordinatesEvaluator(*this); }
	double inline evaluate(eslocal index, size_t timeStep, double temperature, double pressure, double velocity) const
	{
		_values[omp_get_thread_num()][0] = _coordinates[index].x;
		_values[omp_get_thread_num()][1] = _coordinates[index].y;
		_values[omp_get_thread_num()][2] = _coordinates[index].z;
		return _expression[omp_get_thread_num()].evaluate(_values[omp_get_thread_num()]);
	}

	double inline evaluate(const Point &p) const
	{
		_values[omp_get_thread_num()][0] = p.x;
		_values[omp_get_thread_num()][1] = p.y;
		_values[omp_get_thread_num()][2] = p.z;
		return _expression[omp_get_thread_num()].evaluate(_values[omp_get_thread_num()]);
	}

	virtual void store(std::ofstream& os);

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
			Property property = Property::EMPTY);

	TableEvaluator(std::ifstream &is, Property property);

	virtual Evaluator* copy() const { return new TableEvaluator(*this); }
	virtual inline double evaluate(eslocal index, size_t timeStep, double temperature, double pressure, double velocity) const
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

class TableInterpolationEvaluator: public Evaluator {

public:
	TableInterpolationEvaluator(const std::string &name, const std::vector<std::pair<double, double> > &table, Property property = Property::EMPTY);
	TableInterpolationEvaluator(std::ifstream &is, Property property);

	virtual Evaluator* copy() const { return new TableInterpolationEvaluator(*this); }
	virtual inline double evaluate(eslocal index, size_t timeStep, double temperature, double pressure, double velocity) const
	{
		if (temperature < _table[0].first) {
			return _table[0].second;
		}
		for (size_t i = 0; i < _table.size() - 1; i++) {
			if (_table[i].first < temperature && temperature < _table[i + 1].first) {
				double a = _table[i].first  , b = _table[i + 1].first;
				double va = _table[i].second, vb = _table[i + 1].second;
				return va + (vb - va) * (temperature - a) / (b - a);
			}
		}
		return _table.back().second;
	}

	virtual void store(std::ofstream& os);

protected:
	std::vector<std::pair<double, double> > _table;
};

namespace input { class API; }

class ArrayEvaluator: public Evaluator {

	friend class input::API;
public:
	ArrayEvaluator(const std::string &name, std::vector<eslocal> &indices, std::vector<double> &values, eslocal offset, Property property = Property::EMPTY);
	ArrayEvaluator(const std::string &name, eslocal size, eslocal *indices, double *values, eslocal offset, Property property = Property::EMPTY);

	virtual Evaluator* copy() const { return new ArrayEvaluator(*this); }
	virtual inline double evaluate(eslocal index, size_t timeStep, double temperature, double pressure, double velocity) const
	{
		auto it = std::lower_bound(_indices.begin(), _indices.end(), index);
		if (it != _indices.end() && *it == index) {
			return _values[it - _indices.begin()];
		} else {
			ESINFO(ERROR) << "Array evaluator has no specified value for index '" << index << "'";
			return 0;
		}
	}

	virtual void store(std::ofstream& os);

protected:
	std::vector<eslocal> _indices;
	std::vector<double> _values;

};

}



#endif /* SRC_MESH_SETTINGS_EVALUATOR_H_ */
