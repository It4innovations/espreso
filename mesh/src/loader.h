#ifndef LOADING_H_
#define LOADING_H_

#include <algorithm>
#include <string>
#include <sstream>
#include <fstream>
#include <vector>

namespace mesh {

namespace CP {

enum Property {
	DIRICHLET_X,
	DIRICHLET_Y,
	DIRICHLET_Z,
	FORCES_X,
	FORCES_Y,
	FORCES_Z,
	SIZE
};

}

class Loader {

public:
	static size_t getLinesCount(const char* fileName);

private:

	struct TestEOL {
		bool operator()(char c) {
			return c == '\n';
		}
	};
};

class Ansys {

public:
	friend std::ostream& operator<<(std::ostream& os, const Ansys &a);

	Ansys(const char *projectRoot): _projectRoot(projectRoot),
			_elements("ELEMENTS.dat"), _coordinates("COORDINATES.dat"),
			_coordinatesProperty(CP::SIZE) { };

	Ansys(const std::string &projectRoot): _projectRoot(projectRoot),
			_elements("ELEMENTS.dat"), _coordinates("COORDINATES.dat"),
			_coordinatesProperty(CP::SIZE) { };

	std::string coordinates() const
	{
		return _projectRoot + "/" + _coordinates;
	}

	std::string elements() const
	{
		return _projectRoot + "/" + _elements;
	}

	std::string coordinatesProperty(CP::Property property) const
	{
		if (_coordinatesProperty[property].size()) {
			return _projectRoot + "/" + _coordinatesProperty[property];
		} else {
			return "";
		}
	}

	std::string& coordinatesProperty(CP::Property property) {
		return _coordinatesProperty[property];
	}

private:
	std::string _projectRoot;
	std::string _elements;
	std::string _coordinates;
	std::vector<std::string> _coordinatesProperty;
};

}

#endif /* LOADING_H_ */
