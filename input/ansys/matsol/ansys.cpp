
#include "ansys.h"

using namespace esinput;


AnsysMatsol::AnsysMatsol(int argc, char** argv, int rank, int size)
{
	if (argc < 2) {
		if (rank == 0) {
			std::cerr << "Specify the path to an example as the first command line attribute.\n";
		}
		exit(EXIT_FAILURE);
	}
	_path = argv[1];
}

void AnsysMatsol::points(mesh::Coordinates &coordinates)
{
	std::string fileName = _path + "/Model/COORDINATES.dat";

	size_t size = getLinesCount(fileName);
	coordinates.reserve(size);

	std::ifstream file(fileName);
	std::string line;
	mesh::Point point;

	if (file.is_open()) {
		for (size_t c = 0; c < size; c++) {
			getline(file, line, '\n');
			std::stringstream ss(line);

			ss >> point.x >> point.y >> point.z;
			coordinates.add(point, c, c);
		}
		file.close();
	} else {
		std::cerr << "Cannot load mesh from file: " << fileName << "\n";
		exit(EXIT_FAILURE);
	}
}


void AnsysMatsol::elements(std::vector<mesh::Element*> &elements)
{
	int lines;
	std::string settingFile = _path + "/Model/BC/Elasticity/ELEMENT_TYPE.dat";
	std::ifstream setting(settingFile);

	std::string tmp;
	if (setting.is_open()) {
		while (setting) {
			getline(setting, tmp, '\n');
			std::string data = tmp.substr(tmp.find_first_not_of(' '));
			if (data.find("ELEMENT TYPE", 0, 12) == 0) {
				std::stringstream ss(data.substr(data.find("SOLID") + 5, 3));
				eslocal eNumber;
				ss >> eNumber;
				switch (eNumber) {
				case 186: lines = 2; break;
				case 185: lines = 1; break;
				default:
					eslog::error << "Load error: unknown element type\n";
					exit(EXIT_FAILURE);
				}
				break;
			}
		}
	} else {
		std::cerr << "Cannot load element settings from file: " << settingFile << "\n";
		exit(EXIT_FAILURE);
	}

	std::string fileName = _path + "/Model/ELEMENTS.dat";
	elements.resize(getLinesCount(fileName) / lines);

	std::ifstream file(fileName);
	std::vector<std::string> line(lines);

	// 20 is the max of vertices of a element
	// 6 is the parameters number
	eslocal values[20];
	eslocal params[6];
	eslocal value;
	eslocal minIndices = 10000;

	if (file.is_open()) {
		for (eslocal c = 0; c < elements.size(); c++) {
			for (size_t l = 0; l < lines; l++) {
				getline(file, line[l], '\n');
			}

			eslocal n = 0, p = 0;
			for (size_t l = 0; l < line.size(); l++) {
				std::stringstream ss(line[l]);
				while (ss >> value) {
					if (n + p >= 8 && n + p <= 13) {
						params[p++] = value;
					} else {
						values[n++] = value;
					}
				}
			}

			// re-index to zero base
			for (size_t i = 0; i < n; i++) {
				values[i]--;
			}
			elements[c] = AnsysUtils::createElement(values, n);
			elements[c]->setParams(values + n);
		}
		file.close();
	} else {
		std::cerr << "Cannot load mesh from file: " << fileName << "\n";
		exit(EXIT_FAILURE);
	}
}

void AnsysMatsol::boundaryConditions(mesh::Coordinates &coordinates)
{
	std::vector<std::string> conditions = {
		_path + "/Model/BC/Elasticity/NUX.dat",
		_path + "/Model/BC/Elasticity/NUY.dat",
		_path + "/Model/BC/Elasticity/NUZ.dat",
		_path + "/Model/BC/Elasticity/NFX.dat",
		_path + "/Model/BC/Elasticity/NFY.dat",
		_path + "/Model/BC/Elasticity/NFZ.dat"
	};

	for (size_t i = 0; i < coordinates.propertiesSize(); i++) {
		mesh::CoordinatesProperty &property = coordinates.property(i);

		std::ifstream file(conditions[i].c_str());

		if (file.is_open()) {
			eslocal coordinate;
			double value;

			while (file >> coordinate && file.ignore(10, '.') && file >> value) {
				property[coordinate - 1] = value;
			}
			file.close();
		} else {
			std::cout << "Warning: File '" << conditions[i] << "' not found.\n";
		}
	}
}


void AnsysMatsol::clusterBoundaries(mesh::Mesh &mesh, mesh::Boundaries &boundaries)
{
	boundaries.resize(mesh.coordinates().clusterSize());
	for (size_t i = 0; i < mesh.coordinates().clusterSize(); i++) {
		boundaries[i].insert(0);
	}
}

size_t AnsysMatsol::getLinesCount(const std::string &file)
{
	std::ifstream f(file);
	if (f.is_open()) {
		TestEOL test;
		size_t size = count_if(
				std::istreambuf_iterator<char>(f),
				std::istreambuf_iterator<char>(),
				test);

		f.close();
		return size;
	} else {
		std::cerr << "Cannot load file: " << file << "\n";
		exit(EXIT_FAILURE);
	}
}


