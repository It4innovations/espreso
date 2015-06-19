#ifndef LOADING_H_
#define LOADING_H_

#include <algorithm>
#include <string>
#include <sstream>
#include <fstream>
#include <map>

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

	Ansys(const char *projectRoot) :
			_projectRoot(projectRoot), _elements("ELEMENTS.dat"), _coordinates(
					"COORDINATES.dat") {
	}
	std::string coordinates() const{
		return _projectRoot + "/" + _coordinates;
	}

	std::string elements() const{
		return _projectRoot + "/" + _elements;
	}

	const std::map<std::string, std::string>& coordinatesProperties() const{
		return _coordinatesProperties;
	}

	void addCoordinatesProperty(const char* name, const char* file) {
		_coordinatesProperties[name] =_projectRoot+"/"+ file;
	}

private:
	std::string _projectRoot;
	std::string _elements;
	std::string _coordinates;
	std::map<std::string, std::string> _coordinatesProperties;
};

#endif /* LOADING_H_ */
