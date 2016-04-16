
#include "ansys.h"

using namespace espreso::input;

static std::string skip(std::ifstream &is, std::string str)
{
	std::string line;
	while (true) {
		getline(is, line);
		if (line.find(str.c_str(), 0, str.size()) == 0) {
			return line;
		}
	}
}

static eslocal last(std::string &line)
{
	std::stringstream eCode(line.substr(line.find_last_of(',') + 1, line.size()));
	eslocal code;
	eCode >> code;
	return code;
}

void AnsysWorkbench::points(Coordinates &coordinates, size_t &DOFs)
{
	DOFs = 3;
	std::string line;
	Point point;
	size_t id;

	skip(_file, "nblock");
	// skip next line with data format
	getline(_file, line);

	size_t lines = 0;
	while(true) {
		getline(_file, line);
		if (line.find("-1", 0, 2) == 0) {
			break;
		}
		std::stringstream ss(line);
		ss >> id >> point.x >> point.y >> point.z;
		coordinates.add(point, lines, lines);
		lines++;
	}
}


void AnsysWorkbench::elements(std::vector<Element*> &elements)
{
	std::string line;
	size_t lines;

	line = skip(_file, "et");
	switch (last(line)) {
	case 187:
	case 186: lines = 2; break;
	case 185: lines = 1; break;
	default:
		ESINFO(ERROR) << "Load error: unknown element type\n";
	}

	line = skip(_file, "eblock");
	size_t eSize = last(line);

	// (<lineSize>i<numberSize>) e.g. (19i9)
	getline(_file, line);
	std::stringstream format(line);
	eslocal lineSize, numberSize;
	format.ignore(1) >> lineSize;
	format.ignore(1) >> numberSize;

	eslocal values[38], params[11], value, n, p;
	params[0] = 0;
	elements.resize(eSize);
	for (size_t i = 0; i < eSize; i++) {
		p = n = 0;
		for (size_t l = 0; l < lines; l++) {
			char *s = new char[numberSize + 1];
			do {
				_file.get(s, numberSize + 1);
				if (n + p < 11) {
					//params[p++];
					p++;
				} else {
					values[n++] = atol(s);
				}
			} while (_file.peek() != '\n' && _file.peek() != '\r');
			while (_file.peek() == '\n' || _file.peek() == '\r') {
				_file.get();
			}
		}
		for (size_t v = 0; v < n; v++) {
			values[v]--; // re-index
		}
		elements[i] = AnsysUtils::createElement(values, n, params);
	}
}

void AnsysWorkbench::boundaryConditions(Coordinates &coordinates)
{
	std::string line;
	line = skip(_file, "CMBLOCK,SELECTION");

	eslocal value, n = 0, size = last(line);
	getline(_file, line);

	CoordinatesProperty &fx = coordinates.property(FORCES_X);
	while (n < size) {
		getline(_file, line);
		std::stringstream ss(line);
		while (ss >> value) {
			fx[value - 1] = 10;
			n++;
		}
	}

	line = skip(_file, "CMBLOCK,_FIXEDSU");
	n = 0;
	size = last(line);

	getline(_file, line);

	CoordinatesProperty &dx = coordinates.property(DIRICHLET_X);
	CoordinatesProperty &dy = coordinates.property(DIRICHLET_Y);
	CoordinatesProperty &dz = coordinates.property(DIRICHLET_Z);
	while (n < size) {
		getline(_file, line);
		std::stringstream ss(line);
		while (ss >> value) {
			dx[value - 1] = 0;
			dy[value - 1] = 0;
			dz[value - 1] = 0;
			n++;
		}
	}
}


void AnsysWorkbench::clusterBoundaries(Boundaries &boundaries, std::vector<int> &neighbours)
{
	boundaries.resize(mesh.coordinates().clusterSize());
	for (size_t i = 0; i < mesh.coordinates().clusterSize(); i++) {
		boundaries[i].push_back(0);
	}
}


