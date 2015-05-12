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

void Coordinates::computeLocal(size_t part, std::vector<idx_t> &nodeMap, size_t size)
{
	if (_localMappings.size() <= part) {
		_localMappings.resize(part);
	}

	_localMappings[part].clear();
	_localMappings[part].reserve(size);
	for (size_t i = 0; i < nodeMap.size(); i++) {
		if (nodeMap[i] >= 0) {
			_localMappings[part].push_back(i);
		}
	}
}

std::ostream& operator<<(std::ostream& os, const Coordinates &c)
{
	for (size_t i = 0; i < c.size(); i++) {
		os << c._points[i] << "\n";
	}
	return os;
}



