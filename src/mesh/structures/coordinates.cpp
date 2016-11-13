#include "coordinates.h"

namespace espreso {

std::ostream& operator<<(std::ostream& os, const Coordinates &c)
{
	for (size_t i = 0; i < c.clusterSize(); i++)
	{
		os << c._points[i] << "\n";
	}
	return os;
}

}


