#include "coordinates.h"

using namespace mesh;

CoordinatesProperty::CoordinatesProperty(const char* fileName)
{
	std::ifstream file(fileName);
	size_t c = 0;

	if (file.is_open())
	{
		char line[20];
		eslocal coordinate;
		double value;

		while (file >> coordinate && file.ignore(10, '.') && file >> value)
		{
			_mapping[coordinate] = value;
			//std::cout << coordinate << " - " << value << "\n";
		}
		file.close();
	}
	else
	{
		fprintf(stderr, "Cannot load coordinates property from file: %s.\n",
				fileName);
		exit(EXIT_FAILURE);
	}
}

std::ostream& mesh::operator<<(std::ostream& os, const CoordinatesProperty &cp)
{
	std::map<eslocal, double>::const_iterator it = cp._mapping.begin();
	for (; it != cp._mapping.end(); ++it)
	{
		os << it->first << ". " << it->second << "\n";
	}
	return os;
}

Coordinates::Coordinates(const char *fileName) :
		_clusterIndex(1)
{
	_points.resize(Loader::getLinesCount(fileName));
	_clusterIndex[0].resize(_points.size());
	_globalIndex.resize(_points.size());

	std::ifstream file(fileName);
	size_t c = 0;

	if (file.is_open())
	{
		while (c < size() && file >> _points[c])
		{
			_globalIndex[c] = _clusterIndex[0][c] = c++;
		}
		file.close();
	}
	else
	{
		fprintf(stderr, "Cannot load coordinates from file: %s.\n", fileName);
		exit(EXIT_FAILURE);
	}
}

Coordinates::Coordinates(const Coordinates &c) :
		_points(c._points), _clusterIndex(c._clusterIndex), _globalIndex(
				c._globalIndex)
{
	std::map<const std::string, CoordinatesProperty*>::const_iterator it =
			c._coordinatesProperties.begin();
	for (; it != c._coordinatesProperties.end(); ++it)
	{
		_coordinatesProperties[it->first] = new CoordinatesProperty(
				*(it->second));
	}
}

Coordinates & Coordinates::operator=(const Coordinates &c)
{
	if (this != &c)
	{
		_points = c._points;
		_clusterIndex = c._clusterIndex;
		_globalIndex = c._globalIndex;

		std::map<const std::string, CoordinatesProperty*>::const_iterator it =
				_coordinatesProperties.begin();
		for (; it != _coordinatesProperties.end(); ++it)
		{
			delete (it->second);
		}
		_coordinatesProperties.clear();

		it = c._coordinatesProperties.begin();
		for (; it != c._coordinatesProperties.end(); ++it)
		{
			_coordinatesProperties[it->first] = new CoordinatesProperty(
					*(it->second));
		}
	}
	return *this;
}

Coordinates::~Coordinates()
{
	std::map<const std::string, CoordinatesProperty*>::iterator it =
			_coordinatesProperties.begin();
	for (; it != _coordinatesProperties.end(); ++it)
	{
		delete (it->second);
	}
}

void Coordinates::computeLocal(eslocal part, std::vector<eslocal> &nodeMap,
		size_t size)
{
	if (_clusterIndex.size() <= part)
	{
		_clusterIndex.resize(part + 1);
	}

	_clusterIndex[part].clear();
	_clusterIndex[part].reserve(size);
	for (size_t i = 0; i < nodeMap.size(); i++)
	{
		if (nodeMap[i] >= 0)
		{
			_clusterIndex[part].push_back(i);
		}
	}
}

void Coordinates::addCoordinatesProperties(const Ansys &setting)
{
	std::map<std::string, std::string>::const_iterator it =
			setting.coordinatesProperties().begin();
	for (; it != setting.coordinatesProperties().end(); ++it)
	{
		addCoordinatesProperty(it->first,
				new CoordinatesProperty(it->second.c_str()));
	}
}

void Coordinates::printCoordinatesProperties()
{
	std::map<const std::string, CoordinatesProperty*>::const_iterator it =
			_coordinatesProperties.begin();
	for (; it != _coordinatesProperties.end(); ++it)
	{
		std::cout << "Property: " + it->first + "\n";
		std::cout << *(it->second);

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

