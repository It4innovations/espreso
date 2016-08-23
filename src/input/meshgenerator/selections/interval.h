
#ifndef SRC_INPUT_MESHGENERATOR_PARSERS_INTERVAL_H_
#define SRC_INPUT_MESHGENERATOR_PARSERS_INTERVAL_H_

#include <iostream>

#include "esbasis.h"

namespace espreso {

class Interval {

public:

	friend std::ostream& operator<<(std::ostream& os, const Interval& obj);
	friend std::istream& operator>>(std::istream& is, Interval& obj);

	Interval(): _all(true), epsilon(1e-10)
	{
		start[0] = start[1] = start[2] = -epsilon;
		end[0] = end[1] = end[2] = epsilon;
	};

	bool all() const
	{
		return _all;
	}

	double isIn(double x, double y, double z) const
	{
		return all() || (start[0] < x && x < end[0] && start[1] < y && y < end[1] && start[2] < z && z < end[2]);
	}

	double isIn(const Point &p) const
	{
		return isIn(p.x, p.y, p.z);
	}

	double getStart(size_t axis) const
	{
		return start[axis];
	}

	double getEnd(size_t axis) const
	{
		return end[axis];
	}

private:
	bool _all;
	double epsilon;
	double start[3], end[3];
};

}


#endif /* SRC_INPUT_MESHGENERATOR_PARSERS_INTERVAL_H_ */
