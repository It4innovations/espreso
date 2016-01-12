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
	for (size_t i = 0; i < c.clusterSize(); i++)
	{
		os << c._points[i] << "\n";
	}
	return os;
}


