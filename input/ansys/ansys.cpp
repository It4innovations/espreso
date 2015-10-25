
#include "ansys.h"

using namespace esinput;


Ansys::Ansys(int argc, char** argv, int rank, int size)
{
	if (argc < 2) {
		if (rank == 0) {
			std::cerr << "Specify the path to an example as the first command line attribute.\n";
		}
		exit(EXIT_FAILURE);
	}
	_path = argv[1];
}

void Ansys::points(mesh::Coordinates &coordinates)
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


void Ansys::elements(std::vector<mesh::Element*> &elements)
{
	std::string fileName = _path + "/Model/ELEMENTS.dat";
	elements.resize(getLinesCount(fileName));

	std::ifstream file(fileName);
	std::string line;

	// 20 is the max of vertices of a element
	// 6 is the parameters number
	eslocal values[20 + 6], n;
	eslocal value;
	eslocal minIndices = 10000;

	if (file.is_open()) {
		for (eslocal c = 0; c < elements.size(); c++) {
			getline(file, line, '\n');
			std::stringstream ss(line);

			n = 0;
			while (ss >> value) {
				values[n++] = value;
			}
			n -= 6;

			// re-index to zero base
			for (size_t i = 0; i < n; i++) {
				values[i]--;
			}
			elements[c] = createElement(values, n);
			elements[c]->setParams(values + n);
		}
		file.close();
	} else {
		std::cerr << "Cannot load mesh from file: " << fileName << "\n";
		exit(EXIT_FAILURE);
	}
}

void Ansys::boundaryConditions(mesh::Coordinates &coordinates)
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


void Ansys::clusterBoundaries(mesh::Mesh &mesh, mesh::Boundaries &boundaries)
{
	boundaries.resize(mesh.coordinates().size());
	for (size_t i = 0; i < mesh.coordinates().size(); i++) {
		boundaries[i].insert(0);
	}
}


mesh::Element* Ansys::createElement(eslocal *indices, eslocal n)
{
	mesh::Element *e = NULL;
	if (mesh::Tetrahedron4::match(indices, n)) {
		e = new mesh::Tetrahedron4(indices);
	}
	if (mesh::Tetrahedron10::match(indices, n)) {
		e = new mesh::Tetrahedron10(indices);
	}
	if (mesh::Hexahedron8::match(indices, n)) {
		e = new mesh::Hexahedron8(indices);
	}
	if (mesh::Hexahedron20::match(indices, n)) {
		e = new mesh::Hexahedron20(indices);
	}
	if (mesh::Prisma6::match(indices, n)) {
		e = new mesh::Prisma6(indices);
	}
	if (mesh::Prisma15::match(indices, n)) {
		e = new mesh::Prisma15(indices);
	}
	if (mesh::Pyramid5::match(indices, n)) {
		e = new mesh::Pyramid5(indices);
	}
	if (mesh::Pyramid13::match(indices, n)) {
		e = new mesh::Pyramid13(indices);
	}

	if (e == NULL) {
		std::cerr << "Unknown element with indices: ";
		for (eslocal i = 0; i < n; i++) {
			std::cerr << indices[i] << " ";
		}
		std::cerr << "\n";
		exit(EXIT_FAILURE);
	}

	return e;
}

size_t Ansys::getLinesCount(const std::string &file)
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


