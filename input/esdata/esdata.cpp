
#include "esdata.h"

using namespace espreso::input;


void Esdata::points(Coordinates &coordinates, size_t &DOFs)
{
	DOFs = 3; // TODO
	std::stringstream fileName;
	fileName << _path << "/" << _rank << "/coordinates.dat";
	std::ifstream is(fileName.str(), std::ifstream::binary);

	eslocal size;
	esglobal index;
	Point point;

	is.read(reinterpret_cast<char *>(&size), sizeof(eslocal));
	coordinates.reserve(size);
	for (eslocal i = 0; i < size; i++) {
		is.read(reinterpret_cast<char *>(&index), sizeof(esglobal));
		is.read(reinterpret_cast<char *>(&point), Point::size() * sizeof(double));
		coordinates.add(point, i, index);
	}
}


void Esdata::elements(std::vector<Element*> &elements)
{
	std::stringstream fileName;
	fileName << _path << "/" << _rank << "/elements.dat";
	std::ifstream is(fileName.str(), std::ifstream::binary);
	eslocal size, type;
	is.read(reinterpret_cast<char *>(&size), sizeof(eslocal));
	elements.reserve(size);

	for (eslocal i = 0; i < size; i++) {
		is.read(reinterpret_cast<char *>(&type), sizeof(eslocal));
		switch(type) {
		case Tetrahedron4VTKCode:
			elements.push_back(new Tetrahedron4(is));
			break;
		case Tetrahedron10VTKCode:
			elements.push_back(new Tetrahedron10(is));
			break;
		case Pyramid5VTKCode:
			elements.push_back(new Pyramid5(is));
			break;
		case Pyramid13VTKCode:
			elements.push_back(new Pyramid13(is));
			break;
		case Prisma6VTKCode:
			elements.push_back(new Prisma6(is));
			break;
		case Prisma15VTKCode:
			elements.push_back(new Prisma15(is));
			break;
		case Hexahedron8VTKCode:
			elements.push_back(new Hexahedron8(is));
			break;
		case Hexahedron20VTKCode:
			elements.push_back(new Hexahedron20(is));
			break;
		}
	}

	is.close();
}

void Esdata::materials(std::vector<Material> &materials)
{
	// TODO
	materials.resize(1);
}

void Esdata::boundaryConditions(Coordinates &coordinates, std::vector<BoundaryCondition*> &conditions)
{
	std::stringstream fileName;
	fileName << _path << "/" << _rank << "/boundaryConditions.dat";
	std::ifstream is(fileName.str(), std::ifstream::binary);

	eslocal size, cIndex;
	double value;

	for (size_t i = 0; i < coordinates.propertiesSize(); i++) {
		CoordinatesProperty &property = coordinates.property(i);
		is.read(reinterpret_cast<char *>(&size), sizeof(eslocal));
		for (eslocal j = 0; j < size; j++) {
			is.read(reinterpret_cast<char *>(&cIndex), sizeof(eslocal));
			is.read(reinterpret_cast<char *>(&value), sizeof(double));
			property[cIndex] = value;
		}
	}

	size_t counter, DOFs;

	is.read(reinterpret_cast<char *>(&counter), sizeof(size_t));
	for (size_t i = 0; i < counter; i++) {
		is.read(reinterpret_cast<char *>(&DOFs), sizeof(size_t));
		is.read(reinterpret_cast<char *>(&value), sizeof(double));
		conditions.push_back(new NodeCondition(value, ConditionType::DIRICHLET));
		conditions.back()->DOFs().reserve(DOFs);
		for (size_t j = 0; j < DOFs; j++) {
			is.read(reinterpret_cast<char *>(&cIndex), sizeof(eslocal));
			conditions.back()->DOFs().push_back(cIndex);
		}
	}
}


void Esdata::clusterBoundaries(Boundaries &boundaries, std::vector<int> &neighbours)
{
	std::stringstream fileName;
	fileName << _path << "/" << _rank << "/clusterBoundaries.dat";
	std::ifstream is(fileName.str(), std::ifstream::binary);

	std::set<int> neighs;

	boundaries.clear();

	eslocal size, value;

	is.read(reinterpret_cast<char *>(&size), sizeof(eslocal));

	boundaries.resize(size);
	for (size_t i = 0; i < boundaries.size(); i++) {
		is.read(reinterpret_cast<char *>(&size), sizeof(eslocal));
		for (eslocal j = 0; j < size; j++) {
			is.read(reinterpret_cast<char *>(&value), sizeof(eslocal));
			boundaries[i].push_back(value);
			neighs.insert(value);
		}
	}

	neighs.erase(config::env::MPIrank);
	neighbours.insert(neighbours.end(), neighs.begin(), neighs.end());
}

