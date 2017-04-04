
#include <fstream>

#include "../../mesh/elements/line/line2.h"
#include "../../mesh/elements/line/line3.h"

#include "../../mesh/elements/plane/square4.h"
#include "../../mesh/elements/plane/square8.h"
#include "../../mesh/elements/plane/triangle3.h"
#include "../../mesh/elements/plane/triangle6.h"

#include "../../mesh/elements/volume/hexahedron20.h"
#include "../../mesh/elements/volume/hexahedron8.h"
#include "../../mesh/elements/volume/prisma15.h"
#include "../../mesh/elements/volume/prisma6.h"
#include "../../mesh/elements/volume/pyramid13.h"
#include "../../mesh/elements/volume/pyramid5.h"
#include "../../mesh/elements/volume/tetrahedron10.h"
#include "../../mesh/elements/volume/tetrahedron4.h"

#include "../../mesh/structures/mesh.h"
#include "../../mesh/structures/coordinates.h"
#include "../../mesh/structures/region.h"
#include "../../mesh/structures/material.h"
#include "../../mesh/structures/elementtypes.h"
#include "../../mesh/settings/evaluator.h"

#include "espresobinaryformat.h"
#include "../../configuration/environment.h"
#include "../../configuration/input/input.h"

using namespace espreso::input;

void ESPRESOBinaryFormat::load(const ESPRESOInput &configuration, Mesh &mesh, int rank, int size)
{
	auto checkFile = [&] (int cluster, const std::string &file) {
		std::stringstream ss;
		ss << configuration.path << "/" << cluster << "/" << file;
		std::ifstream is (ss.str());
		if (!is.good()) {
			ESINFO(GLOBAL_ERROR) << "Old ESPRESO binary format. File '" << ss.str() << "' is missing.";
		}
	};

	ESINFO(OVERVIEW) << "Load mesh from ESPRESO binary format from directory " << configuration.path;
	ESPRESOBinaryFormat esdata(configuration, mesh, rank, size);
	std::ifstream is(configuration.path + "/description.txt");
	if (!is.good()) {
		ESINFO(GLOBAL_ERROR) << "Old ESPRESO binary format. File '" << configuration.path << "/description.txt' is missing.";
	}
	int clusters;
	is >> clusters;
	if (clusters != environment->MPIsize) {
		ESINFO(GLOBAL_ERROR) << "Incorrect number of MPI processes (" << environment->MPIsize << "). Should be " << clusters;
	}



	for (int cluster = 0; cluster < clusters; cluster++) {
		checkFile(cluster, "coordinates.dat");
		checkFile(cluster, "elements.dat");
		checkFile(cluster, "materials.dat");
		checkFile(cluster, "regions.dat");
		checkFile(cluster, "boundaries.dat");
	}

	esdata.fill();
}

void ESPRESOBinaryFormat::points(Coordinates &coordinates)
{
	std::stringstream fileName;
	fileName << _esdata.path << "/" << _rank << "/coordinates.dat";
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

void ESPRESOBinaryFormat::elements(std::vector<size_t> &bodies, std::vector<Element*> &elements, std::vector<Element*> &faces, std::vector<Element*> &edges)
{
	std::stringstream fileName;
	fileName << _esdata.path << "/" << _rank << "/elements.dat";
	std::ifstream is(fileName.str(), std::ifstream::binary);
	eslocal size;

	is.read(reinterpret_cast<char *>(&size), sizeof(eslocal));
	elements.reserve(size);
	addElements(is, elements, size);

	is.close();
	bodies = { 0, elements.size() };
}

void ESPRESOBinaryFormat::materials(std::vector<Material*> &materials)
{
	std::stringstream fileName;
	fileName << _esdata.path << "/" << _rank << "/materials.dat";
	std::ifstream is(fileName.str(), std::ifstream::binary);

	eslocal size;
	is.read(reinterpret_cast<char *>(&size), sizeof(eslocal));
	for (eslocal i = 0; i < size; i++) {
		materials.push_back(new Material(mesh.coordinates()));
		materials.back()->load(is);
	}
	is.close();
}

void ESPRESOBinaryFormat::regions(
		std::vector<Evaluator*> &evaluators,
		std::vector<Region*> &regions,
		std::vector<Element*> &elements,
		std::vector<Element*> &faces,
		std::vector<Element*> &edges,
		std::vector<Element*> &nodes)
{
	std::stringstream fileName;
	fileName << _esdata.path << "/" << _rank << "/regions.dat";
	std::ifstream is(fileName.str(), std::ifstream::binary);

	eslocal size;
	int index;

	is.read(reinterpret_cast<char *>(&size), sizeof(eslocal));
	faces.reserve(size);
	addElements(is, faces, size);

	is.read(reinterpret_cast<char *>(&size), sizeof(eslocal));
	edges.reserve(size);
	addElements(is, edges, size);

	is.read(reinterpret_cast<char *>(&size), sizeof(eslocal));
	for (eslocal i = 0; i < size; i++) {
		evaluators.push_back(Evaluator::create(is, mesh.coordinates()));
	}

	is.read(reinterpret_cast<char *>(&size), sizeof(eslocal));
	for (eslocal i = 0; i < size; i++) {
		if (i > 1) {
			regions.push_back(new Region(ElementType::NODES));
			eslocal length;
			is.read(reinterpret_cast<char *>(&length), sizeof(eslocal));
			char *buffer = new char[length];
			is.read(buffer, length);
			regions.back()->name = std::string(buffer, buffer + length);
			delete[] buffer;
			is.read(reinterpret_cast<char *>(&regions.back()->eType), sizeof(ElementType));
		}

		eslocal steps;
		is.read(reinterpret_cast<char *>(&steps), sizeof(eslocal));
		regions[i]->settings.resize(steps);
		for (eslocal step = 0; step < steps; step++) {
			eslocal properties;
			is.read(reinterpret_cast<char *>(&properties), sizeof(eslocal));
			for (eslocal property = 0; property < properties; property++) {
				int p, index;
				eslocal eSize;
				is.read(reinterpret_cast<char *>(&p), sizeof(int));
				is.read(reinterpret_cast<char *>(&eSize), sizeof(eslocal));
				for (eslocal j = 0; j < eSize; j++) {
					is.read(reinterpret_cast<char *>(&index), sizeof(int));
					regions[i]->settings[step][(Property)p].push_back(evaluators[index]);
				}
			}
		}
	}

	is.read(reinterpret_cast<char *>(&size), sizeof(eslocal));
	ESTEST(MANDATORY) << "Invalid size element regions" << (size == (int)elements.size() ? TEST_PASSED : TEST_FAILED);
	for (size_t i = 0; i < elements.size(); i++) {
		is.read(reinterpret_cast<char *>(&size), sizeof(eslocal));
		for (eslocal r = 0; r < size; r++) {
			is.read(reinterpret_cast<char *>(&index), sizeof(int));
			regions[index]->elements().push_back(elements[i]);
		}
	}

	is.read(reinterpret_cast<char *>(&size), sizeof(eslocal));
	ESTEST(MANDATORY) << "Invalid size of faces regions" << (size == (int)faces.size() ? TEST_PASSED : TEST_FAILED);
	for (size_t i = 0; i < faces.size(); i++) {
		is.read(reinterpret_cast<char *>(&size), sizeof(eslocal));
		for (eslocal r = 0; r < size; r++) {
			is.read(reinterpret_cast<char *>(&index), sizeof(int));
			regions[index]->elements().push_back(faces[i]);
		}
	}

	is.read(reinterpret_cast<char *>(&size), sizeof(eslocal));
	ESTEST(MANDATORY) << "Invalid size of edge regions" << (size == (int)edges.size() ? TEST_PASSED : TEST_FAILED);
	for (size_t i = 0; i < edges.size(); i++) {
		is.read(reinterpret_cast<char *>(&size), sizeof(eslocal));
		for (eslocal r = 0; r < size; r++) {
			is.read(reinterpret_cast<char *>(&index), sizeof(int));
			regions[index]->elements().push_back(edges[i]);
		}
	}

	is.read(reinterpret_cast<char *>(&size), sizeof(eslocal));
	ESTEST(MANDATORY) << "Invalid size of node regions" << (size == (int)nodes.size() ? TEST_PASSED : TEST_FAILED);
	for (size_t i = 0; i < nodes.size(); i++) {
		is.read(reinterpret_cast<char *>(&size), sizeof(eslocal));
		for (eslocal r = 0; r < size; r++) {
			is.read(reinterpret_cast<char *>(&index), sizeof(int));
			regions[index]->elements().push_back(nodes[i]);
		}
	}

	is.close();
}

bool ESPRESOBinaryFormat::partitiate(const std::vector<Element*> &nodes, std::vector<eslocal> &partsPtrs, std::vector<std::vector<Element*> > &fixPoints, std::vector<Element*> &corners)
{
	mesh.partitiate(_esdata.domains);
	return true;
}

void ESPRESOBinaryFormat::neighbours(std::vector<Element*> &nodes, std::vector<int> &neighbours, const std::vector<Element*> &faces, const std::vector<Element*> &edges)
{
	std::stringstream fileName;
	fileName << _esdata.path << "/" << _rank << "/boundaries.dat";
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

	neighs.erase(environment->MPIrank);
	neighbours.insert(neighbours.end(), neighs.begin(), neighs.end());
}

