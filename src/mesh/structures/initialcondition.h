
#ifndef SRC_MESH_STRUCTURES_INITIALCONDITION_H_
#define SRC_MESH_STRUCTURES_INITIALCONDITION_H_

#include "esbasis.h"

namespace espreso {

class InitialCondition {

public:
	virtual InitialCondition* copy() =0;

	virtual double get(const Point &p) =0;
	virtual ~InitialCondition() {};
};

class IntervalInitialization: public InitialCondition {

public:
	double get(const Point &p)
	{
		for (size_t i = 0; i < intervals.size(); i++) {
			if (intervals[i].first.isIn(p)) {
				return intervals[i].second.evaluate(p.x, p.y, p.z);
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
		intervals.push_back(std::make_pair(interval, Expression(expression)));
	}

protected:
	std::vector<std::pair<Interval, Expression> > intervals;
};

}



#endif /* SRC_MESH_STRUCTURES_INITIALCONDITION_H_ */
