
#ifndef SRC_MESH_STRUCTURES_INITIALCONDITION_H_
#define SRC_MESH_STRUCTURES_INITIALCONDITION_H_

#include "esbasis.h"

namespace espreso {

enum class InitialConditionType {
	TRANSLATION_MOTION_X,
	TRANSLATION_MOTION_Y,
	TRANSLATION_MOTION_Z,
	HEAT_SOURCE,
	INITIAL_CONDITION_SIZE
};

class InitialCondition {

public:
	static int index(InitialConditionType type)
	{
		return static_cast<int>(type);
	}

	virtual InitialCondition* copy() =0;

	virtual double get(const Point &p) const =0;
	virtual ~InitialCondition() {};
};

class ExpressionInitialization: public InitialCondition {

public:
	double get(const Point &p) const
	{
		return _expression.evaluate(p.x, p.y, p.z);
	}

	ExpressionInitialization* copy()
	{
		return new ExpressionInitialization(*this);
	}

	ExpressionInitialization(const std::string &expression): _expression(expression) {};

protected:
	Expression _expression;
};

class IntervalInitialization: public InitialCondition {

public:
	double get(const Point &p) const
	{
		for (size_t i = 0; i < _intervals.size(); i++) {
			if (_intervals[i].first.isIn(p)) {
				return _intervals[i].second.evaluate(p.x, p.y, p.z);
			}
		}
		ESINFO(ERROR) << "Not filled initial condition for point " << p;
		return 0;
	}

	IntervalInitialization* copy()
	{
		return new IntervalInitialization(*this);
	}

	void add(const Interval &interval, const std::string &expression)
	{
		_intervals.push_back(std::make_pair(interval, Expression(expression)));
	}

protected:
	std::vector<std::pair<Interval, Expression> > _intervals;
};

}



#endif /* SRC_MESH_STRUCTURES_INITIALCONDITION_H_ */
