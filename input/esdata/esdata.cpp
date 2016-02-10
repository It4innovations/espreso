
#include "esdata.h"

using namespace esinput;

Esdata::Esdata(const Options &options, int rank, int size): _rank(rank), _size(size)
{
	_path = options.path;
}


void Esdata::points(mesh::Coordinates &coordinates)
{
	std::stringstream fileName;
	fileName << _path << "/" << _rank << "/coordinates.dat";
	std::ifstream is(fileName.str(), std::ifstream::binary);

	eslocal size;
	esglobal index;
	mesh::Point point;

	is.read(reinterpret_cast<char *>(&size), sizeof(eslocal));
	coordinates.reserve(size);
	for (eslocal i = 0; i < size; i++) {
		is.read(reinterpret_cast<char *>(&index), sizeof(esglobal));
		is.read(reinterpret_cast<char *>(&point), mesh::Point::size() * sizeof(double));
		coordinates.add(point, i, index);
	}
}


void Esdata::elements(std::vector<mesh::Element*> &elements)
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
			elements.push_back(new mesh::Tetrahedron4(is));
			break;
		case Tetrahedron10VTKCode:
			elements.push_back(new mesh::Tetrahedron10(is));
			break;
		case Pyramid5VTKCode:
			elements.push_back(new mesh::Pyramid5(is));
			break;
		case Pyramid13VTKCode:
			elements.push_back(new mesh::Pyramid13(is));
			break;
		case Prisma6VTKCode:
			elements.push_back(new mesh::Prisma6(is));
			break;
		case Prisma15VTKCode:
			elements.push_back(new mesh::Prisma15(is));
			break;
		case Hexahedron8VTKCode:
			elements.push_back(new mesh::Hexahedron8(is));
			break;
		case Hexahedron20VTKCode:
			elements.push_back(new mesh::Hexahedron20(is));
			break;
		}
	}

	is.close();
}

void Esdata::boundaryConditions(mesh::Coordinates &coordinates)
{
	std::stringstream fileName;
	fileName << _path << "/" << _rank << "/boundaryConditions.dat";
	std::ifstream is(fileName.str(), std::ifstream::binary);

	eslocal size, cIndex;
	double value;

	for (size_t i = 0; i < coordinates.propertiesSize(); i++) {
		mesh::CoordinatesProperty &property = coordinates.property(i);
		is.read(reinterpret_cast<char *>(&size), sizeof(eslocal));
		for (eslocal j = 0; j < size; j++) {
			is.read(reinterpret_cast<char *>(&cIndex), sizeof(eslocal));
			is.read(reinterpret_cast<char *>(&value), sizeof(double));
			property[cIndex] = value;
		}
	}
}


void Esdata::clusterBoundaries(mesh::Mesh &mesh, mesh::Boundaries &boundaries)
{
	std::stringstream fileName;
	fileName << _path << "/" << _rank << "/clusterBoundaries.dat";
	std::ifstream is(fileName.str(), std::ifstream::binary);

	boundaries.clear();

	eslocal size, value;

	is.read(reinterpret_cast<char *>(&size), sizeof(eslocal));

	boundaries.resize(size);
	for (size_t i = 0; i < boundaries.size(); i++) {
		is.read(reinterpret_cast<char *>(&size), sizeof(eslocal));
		for (eslocal j = 0; j < size; j++) {
			is.read(reinterpret_cast<char *>(&value), sizeof(eslocal));
			boundaries[i].insert(value);
		}
	}
}

