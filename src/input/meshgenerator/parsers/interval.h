
#ifndef SRC_INPUT_MESHGENERATOR_PARSERS_INTERVAL_H_
#define SRC_INPUT_MESHGENERATOR_PARSERS_INTERVAL_H_

#include "esbasis.h"
#include "esmesh.h"

namespace espreso {
namespace input {

class Interval {

public:

	friend std::ostream& operator<<(std::ostream& os, const Interval& obj);
	friend std::istream& operator>>(std::istream& is, Interval& obj);

	Interval()
	: start{0, 0, 0}, end{0, 0, 0}, excludeStart{false, false, false}, excludeEnd{false, false, false} {};


	double isIn(double x, double y, double z) const
	{
		return
			(((start[0] < x) || (!excludeStart[0] && start[0] == x)) && ((x < end[0]) || (!excludeEnd[0] && x == end[0]))) &&
			(((start[1] < y) || (!excludeStart[1] && start[1] == y)) && ((y < end[1]) || (!excludeEnd[1] && y == end[1]))) &&
			(((start[2] < z) || (!excludeStart[2] && start[2] == z)) && ((z < end[2]) || (!excludeEnd[2] && z == end[2])));
	}

	double isIn(const Point3D &p) const
	{
		return isIn(p.x, p.y, p.z);
	}

private:
	double start[3], end[3];
	bool excludeStart[3], excludeEnd[3];
};


}
}


#endif /* SRC_INPUT_MESHGENERATOR_PARSERS_INTERVAL_H_ */
