#include "coordinates.h"

Coordinates::Coordinates(const char *fileName): _offset(1)
{
	_points.resize(Loader::getLinesCount(fileName));

	std::ifstream file(fileName);
	size_t c = 0;

	if (file.is_open()) {
		while(c < size() && file >> _points[c++]);
		file.close();
	} else {
		fprintf(stderr, "Cannot load coordinates from file: %s.\n", fileName);
		exit(EXIT_FAILURE);
	}
}

std::ostream& operator<<(std::ostream& os, const Coordinates &c)
{
	for (size_t i = 0; i < c.size(); i++) {
		os << c._points[i] << "\n";
	}
	return os;
}



