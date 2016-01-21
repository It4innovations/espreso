
#include "ansys.h"

using namespace esinput;


AnsysWorkbench::AnsysWorkbench(int argc, char** argv, size_t index, size_t size)
{
	if (argc < 2) {
		if (index == 0) {
			std::cerr << "Specify the path to an example as the first command line attribute.\n";
		}
		exit(EXIT_FAILURE);
	}
	_path = argv[1];
}

//static size_t getCoordinateSize(std::ifstream &is)
//{
//	std::streampos start = is.tellg();
//	std::string line;
//	size_t lineCount = 0;
//	while(true) {
//		getline(is, line);
//		if (line.find("-1", 0, 2) != std::string::npos) {
//			break;
//		}
//		lineCount++;
//	}
//	is.seekg(start);
//	return lineCount;
//}

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

void AnsysWorkbench::points(mesh::Coordinates &coordinates)
{
	std::string line;
	mesh::Point point;
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


void AnsysWorkbench::elements(std::vector<mesh::Element*> &elements)
{
	std::string line;
	size_t lines;

	line = skip(_file, "et");
	switch (last(line)) {
	case 187:
	case 186: lines = 2; break;
	case 185: lines = 1; break;
	default:
		eslog::error << "Load error: unknown element type\n";
		exit(EXIT_FAILURE);
	}

	line = skip(_file, "eblock");
	size_t eSize = last(line);
	getline(_file, line);

	eslocal values[38], value, n, p;
	elements.resize(eSize);
	for (size_t i = 0; i < eSize; i++) {
		p = n = 0;
		for (size_t l = 0; l < lines; l++) {
			getline(_file, line);
			std::stringstream ss(line);
			while (ss >> value) {
				if (n + p < 11) {
					p++; // TODO: do not skip parameters
				} else {
					values[n++] = value;
				}
			}
		}
		for (size_t v = 0; v < n; v++) {
			values[v]--; // re-index
		}
		elements[i] = AnsysUtils::createElement(values, n);
	}
}

void AnsysWorkbench::boundaryConditions(mesh::Coordinates &coordinates)
{
	std::string line;
	line = skip(_file, "CMBLOCK");
	eslocal value, n = 0, size = last(line);
	getline(_file, line);

	mesh::CoordinatesProperty &dx = coordinates.property(mesh::DIRICHLET_X);
	mesh::CoordinatesProperty &dy = coordinates.property(mesh::DIRICHLET_Y);
	mesh::CoordinatesProperty &dz = coordinates.property(mesh::DIRICHLET_Z);
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


void AnsysWorkbench::clusterBoundaries(mesh::Mesh &mesh, mesh::Boundaries &boundaries)
{
	boundaries.resize(mesh.coordinates().clusterSize());
	for (size_t i = 0; i < mesh.coordinates().clusterSize(); i++) {
		boundaries[i].insert(0);
	}
}


