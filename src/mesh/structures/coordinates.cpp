#include "coordinates.h"

using namespace espreso;

std::ostream& espreso::operator<<(std::ostream& os, const Coordinates &c)
{
	for (size_t i = 0; i < c.clusterSize(); i++)
	{
		os << c._points[i] << "\n";
	}
	return os;
}


