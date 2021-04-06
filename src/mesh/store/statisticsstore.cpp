
#include "statisticsstore.h"

#include <limits>

using namespace espreso;

Statistics::Statistics()
{
	reset();
}

void Statistics::reset()
{
	min = std::numeric_limits<double>::max();
	max = -std::numeric_limits<double>::max();
	avg = norm = 0;
	absmin = std::numeric_limits<double>::max();
	absmax = 0;
}



