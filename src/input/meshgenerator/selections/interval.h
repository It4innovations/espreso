
#ifndef SRC_INPUT_MESHGENERATOR_PARSERS_INTERVAL_H_
#define SRC_INPUT_MESHGENERATOR_PARSERS_INTERVAL_H_

#include <iostream>

#include "esbasis.h"

namespace espreso {

class Interval {

public:

	friend std::ostream& operator<<(std::ostream& os, const Interval& obj);
	friend std::istream& operator>>(std::istream& is, Interval& obj);

	Interval()
	: _all(true),
	  start{0, 0, 0}, end{0, 0, 0},
	  excludeStart{false, false, false}, excludeEnd{false, false, false} {}

	Interval(double sx, double ex, double sy, double ey, double sz, double ez)
	: _all(false),
	  start{sx, sy, sz}, end{ex, ey, ez},
	  excludeStart{false, false, false}, excludeEnd{false, false, false} {}

	bool all() const
	{
		return _all;
	}

	double start[3], end[3];
	bool excludeStart[3], excludeEnd[3];

private:
	bool _all;

};

}


#endif /* SRC_INPUT_MESHGENERATOR_PARSERS_INTERVAL_H_ */
