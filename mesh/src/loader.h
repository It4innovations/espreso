#ifndef LOADING_H_
#define LOADING_H_

#include <algorithm>
#include <string>
#include <sstream>
#include <fstream>

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


private:
	std::string _projectRoot;
	std::string _elements;
	std::string _coordinates;
};

#endif /* LOADING_H_ */
