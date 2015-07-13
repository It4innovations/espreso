#include "coordinates.h"

using namespace mesh;

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
		std::cout << "File '" << fileName << "' not found.\n";
	}
}

std::ostream& mesh::operator<<(std::ostream& os, const CoordinatesProperty &cp)
{
	std::map<eslocal, double>::const_iterator it;
	for (it = cp._mapping.begin(); it != cp._mapping.end(); ++it) {
		os << it->first << ". " << it->second << "\n";
	}
	return os;
}

void Coordinates::readFromFile(const char *fileName)
{
	_points.resize(Loader::getLinesCount(fileName));
	_clusterIndex[0].resize(_points.size());
	_globalIndex.resize(_points.size());

	std::ifstream file(fileName);
	size_t c = 0;

	if (file.is_open()) {
		while (c < size() && file >> _points[c]) {
			_globalIndex[c] = _clusterIndex[0][c] = c++;
		}
		file.close();
	} else {
		fprintf(stderr, "Cannot load coordinates from file: %s.\n", fileName);
		exit(EXIT_FAILURE);
	}
}

Coordinates::Coordinates(const char *fileName): _clusterIndex(1), _property(CP::SIZE)
{
	readFromFile(fileName);
}

Coordinates::Coordinates(const Ansys &a): _clusterIndex(1), _property(CP::SIZE)
{
	readFromFile(a.coordinates().c_str());
	_property[CP::DIRICHLET_X].load(a.coordinatesProperty(CP::DIRICHLET_X).c_str());
	_property[CP::DIRICHLET_Y].load(a.coordinatesProperty(CP::DIRICHLET_Y).c_str());
	_property[CP::DIRICHLET_Z].load(a.coordinatesProperty(CP::DIRICHLET_Z).c_str());
	_property[CP::FORCES_X].load(a.coordinatesProperty(CP::FORCES_X).c_str());
	_property[CP::FORCES_Y].load(a.coordinatesProperty(CP::FORCES_Y).c_str());
	_property[CP::FORCES_Z].load(a.coordinatesProperty(CP::FORCES_Z).c_str());
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
	for (size_t i = 0; i < c.size(); i++)
	{
		os << c._points[i] << "\n";
	}
	return os;
}


