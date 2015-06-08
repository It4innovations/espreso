#include "coordinates.h"

using namespace mesh;

Coordinates::Coordinates(const char *fileName): _clusterIndex(1)
{
	_points.resize(Loader::getLinesCount(fileName));
	_clusterIndex[0].resize(_points.size());
	_globalIndex.resize(_points.size());

	std::ifstream file(fileName);
	size_t c = 0;

	if (file.is_open()) {
		while(c < size() && file >> _points[c])
		{
			_globalIndex[c] = _clusterIndex[0][c] = c++;
		}
		file.close();
	} else {
		fprintf(stderr, "Cannot load coordinates from file: %s.\n", fileName);
		exit(EXIT_FAILURE);
	}
}

void Coordinates::computeLocal(eslocal part, std::vector<eslocal> &nodeMap, size_t size)
{
	if (_clusterIndex.size() <= part) {
		_clusterIndex.resize(part + 1);
	}

	_clusterIndex[part].clear();
	_clusterIndex[part].reserve(size);
	for (size_t i = 0; i < nodeMap.size(); i++) {
		if (nodeMap[i] >= 0) {
			_clusterIndex[part].push_back(i);
		}
	}
}

std::ostream& mesh::operator<<(std::ostream& os, const Coordinates &c)
{
	for (size_t i = 0; i < c.size(); i++) {
		os << c._points[i] << "\n";
	}
	return os;
}

