
#include "esdata.h"

using namespace espreso::input;


void Esdata::points(Coordinates &coordinates)
{
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

static void addElements(std::ifstream &is, std::vector<espreso::Element*> &elements, size_t number)
{
	eslocal type;
	for (size_t i = 0; i < number; i++) {
		is.read(reinterpret_cast<char *>(&type), sizeof(eslocal));
		switch(type) {
		case Tetrahedron4VTKCode:
			elements.push_back(new espreso::Tetrahedron4(is));
			break;
		case Tetrahedron10VTKCode:
			elements.push_back(new espreso::Tetrahedron10(is));
			break;
		case Pyramid5VTKCode:
			elements.push_back(new espreso::Pyramid5(is));
			break;
		case Pyramid13VTKCode:
			elements.push_back(new espreso::Pyramid13(is));
			break;
		case Prisma6VTKCode:
			elements.push_back(new espreso::Prisma6(is));
			break;
		case Prisma15VTKCode:
			elements.push_back(new espreso::Prisma15(is));
			break;
		case Hexahedron8VTKCode:
			elements.push_back(new espreso::Hexahedron8(is));
			break;
		case Hexahedron20VTKCode:
			elements.push_back(new espreso::Hexahedron20(is));
			break;
		case Square4VTKCode:
			elements.push_back(new espreso::Square4(is));
			break;
		case Square8VTKCode:
			elements.push_back(new espreso::Square8(is));
			break;
		case Triangle3VTKCode:
			elements.push_back(new espreso::Triangle3(is));
			break;
		case Triangle6VTKCode:
			elements.push_back(new espreso::Triangle6(is));
			break;
		case Line2VTKCode:
			elements.push_back(new espreso::Line2(is));
			break;
		case Line3VTKCode:
			elements.push_back(new espreso::Line3(is));
			break;
		}
	}
};

void Esdata::elements(std::vector<Element*> &elements, std::vector<Element*> &faces, std::vector<Element*> &edges)
{
	std::stringstream fileName;
	fileName << _path << "/" << _rank << "/elements.dat";
	std::ifstream is(fileName.str(), std::ifstream::binary);
	eslocal size;

	is.read(reinterpret_cast<char *>(&size), sizeof(eslocal));
	elements.reserve(size);
	addElements(is, elements, size);

	is.close();
}

void Esdata::materials(std::vector<Material> &materials)
{
	std::stringstream fileName;
	fileName << _path << "/" << _rank << "/materials.dat";
	std::ifstream is(fileName.str(), std::ifstream::binary);

	eslocal size;
	is.read(reinterpret_cast<char *>(&size), sizeof(eslocal));
	for (eslocal i = 0; i < size; i++) {
		materials.push_back(Material(is, mesh.coordinates()));
	}
	is.close();
}

void Esdata::regions(
		std::vector<Evaluator*> &evaluators,
		std::vector<Region> &regions,
		std::vector<Element*> &elements,
		std::vector<Element*> &faces,
		std::vector<Element*> &edges,
		std::vector<Element*> &nodes)
{
	std::stringstream fileName;
	fileName << _path << "/" << _rank << "/settings.dat";
	std::ifstream is(fileName.str(), std::ifstream::binary);

	eslocal size;

	is.read(reinterpret_cast<char *>(&size), sizeof(eslocal));
	faces.reserve(size);
	addElements(is, elements, size);

	is.read(reinterpret_cast<char *>(&size), sizeof(eslocal));
	edges.reserve(size);
	addElements(is, elements, size);

	is.read(reinterpret_cast<char *>(&size), sizeof(eslocal));
	for (eslocal i = 0; i < size; i++) {
		evaluators.push_back(Evaluator::create(is, mesh.coordinates()));
	}

	is.read(reinterpret_cast<char *>(&size), sizeof(eslocal));
	ESTEST(MANDATORY) << "Invalid size of element settings" << (size == (int)elements.size() ? TEST_PASSED : TEST_FAILED);
	for (eslocal i = 0; i < size; i++) {
		elements[i]->settings().load(is, evaluators);
	}

	is.read(reinterpret_cast<char *>(&size), sizeof(eslocal));
	ESTEST(MANDATORY) << "Invalid size of faces settings" << (size == (int)faces.size() ? TEST_PASSED : TEST_FAILED);
	for (eslocal i = 0; i < size; i++) {
		faces[i]->settings().load(is, evaluators);
	}

	is.read(reinterpret_cast<char *>(&size), sizeof(eslocal));
	ESTEST(MANDATORY) << "Invalid size of edge settings" << (size == (int)edges.size() ? TEST_PASSED : TEST_FAILED);
	for (eslocal i = 0; i < size; i++) {
		edges[i]->settings().load(is, evaluators);
	}

	is.read(reinterpret_cast<char *>(&size), sizeof(eslocal));
	ESTEST(MANDATORY) << "Invalid size of node settings" << (size == (int)nodes.size() ? TEST_PASSED : TEST_FAILED);
	for (eslocal i = 0; i < size; i++) {
		nodes[i]->settings().load(is, evaluators);
	}

	is.close();
}


void Esdata::neighbours(std::vector<Element*> &nodes, std::vector<int> &neighbours)
{
	std::stringstream fileName;
	fileName << _path << "/" << _rank << "/boundaries.dat";
	std::ifstream is(fileName.str(), std::ifstream::binary);

	std::set<int> neighs;

	eslocal size, nSize, value;

	is.read(reinterpret_cast<char *>(&nSize), sizeof(eslocal));
	ESTEST(MANDATORY) << "Invalid node size of cluster boundaries" << (nSize == (int)nodes.size() ? TEST_PASSED : TEST_FAILED);
	for (eslocal i = 0; i < nSize; i++) {
		is.read(reinterpret_cast<char *>(&size), sizeof(eslocal));
		for (eslocal j = 0; j < size; j++) {
			is.read(reinterpret_cast<char *>(&value), sizeof(eslocal));
			nodes[i]->clusters().push_back(value);
			neighs.insert(value);
		}
	}

	neighs.erase(config::env::MPIrank);
	neighbours.insert(neighbours.end(), neighs.begin(), neighs.end());
}

