#include "coordinates.h"

using namespace espreso;

void CoordinatesProperty::load(const char* fileName)
{
	std::ifstream file(fileName);

	if (file.is_open()) {
		eslocal coordinate;
		double value;

		while (file >> coordinate && file.ignore(10, '.') && file >> value) {
			_mapping[coordinate - 1] = value;
		}
		file.close();
	} else {
		ESINFO(ERROR) << "File '" << fileName << "' not found";
	}
}

std::ostream& espreso::operator<<(std::ostream& os, const CoordinatesProperty &cp)
{
	std::map<eslocal, double>::const_iterator it;
	for (it = cp._mapping.begin(); it != cp._mapping.end(); ++it) {
		os << it->first << ". " << it->second << "\n";
	}
	return os;
}

std::ostream& espreso::operator<<(std::ostream& os, const Coordinates &c)
{
	for (size_t i = 0; i < c.clusterSize(); i++)
	{
		os << c._points[i] << "\n";
	}
	return os;
}


