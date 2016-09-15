
#ifndef SRC_MESH_SETTINGS_EVALUATOR_H_
#define SRC_MESH_SETTINGS_EVALUATOR_H_

#include "esbasis.h"
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
		COORDINATE
	};

public:
	static Evaluator* create(std::ifstream &is, const Coordinates &coordinates);

	virtual Evaluator* copy() const { return new Evaluator(*this); }
	virtual double evaluate(size_t index) const { return 0; }
	virtual ~Evaluator() {};

	virtual void store(std::ofstream& os)
	{
		Type type = Type::DEFAULT;
		os.write(reinterpret_cast<const char *>(&type), sizeof(Evaluator::Type));
	}
};

class ConstEvaluator: public Evaluator {

public:
	ConstEvaluator(double value): _value(value) {};
	ConstEvaluator(std::ifstream &is) { is.read(reinterpret_cast<char *>(&_value), sizeof(double)); }

	virtual Evaluator* copy() const { return new ConstEvaluator(*this); }
	double evaluate(size_t index) const { return _value; }

	virtual void store(std::ofstream& os)
	{
		Type type = Type::CONST;
		os.write(reinterpret_cast<const char *>(&type), sizeof(Evaluator::Type));
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
	ExpressionEvaluator(const std::string &expression, std::vector<std::string> variables)
	: _expression(config::env::CILK_NWORKERS, Expression(expression, variables)),
	  _values(config::env::CILK_NWORKERS, std::vector<double>(variables.size())) {};

	ExpressionEvaluator(std::ifstream &is, std::vector<std::string> variables)
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
	virtual double evaluate(size_t index) const =0;

	std::vector<Expression> _expression;
	mutable std::vector<std::vector<double> >_values;
};

class CoordinatesEvaluator: public ExpressionEvaluator {

public:
	CoordinatesEvaluator(const std::string &expression, const Coordinates &coordinates)
	: ExpressionEvaluator(expression, { "x", "y", "z" }), _coordinates(coordinates) {};

	CoordinatesEvaluator(std::ifstream &is, const Coordinates &coordinates)
	: ExpressionEvaluator(is, { "x", "y", "z" }), _coordinates(coordinates) {};

	virtual Evaluator* copy() const { return new CoordinatesEvaluator(*this); }
	double evaluate(size_t index) const
	{
		_values[__cilkrts_get_worker_number()][0] = _coordinates[index].x;
		_values[__cilkrts_get_worker_number()][1] = _coordinates[index].y;
		_values[__cilkrts_get_worker_number()][2] = _coordinates[index].z;
		return _expression[__cilkrts_get_worker_number()].evaluate(_values[__cilkrts_get_worker_number()]);
	}

	virtual void store(std::ofstream& os)
	{
		Type type = Type::COORDINATE;
		os.write(reinterpret_cast<const char *>(&type), sizeof(Evaluator::Type));
		os.write(_expression[0].expression().c_str(), _expression[0].expression().size());
	}

private:
	const Coordinates &_coordinates;
};

}



#endif /* SRC_MESH_SETTINGS_EVALUATOR_H_ */
