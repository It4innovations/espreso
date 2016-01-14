
#include "generator.h"
#include "factory/factory.h"

using namespace esinput;


MeshGenerator::MeshGenerator(int argc, char** argv, size_t index, size_t size)
{
	_generator = MeshFactory::create(argc, argv, index, size);
}

bool MeshGenerator::manualPartition()
{
	return _generator->manualPartition();
}

void MeshGenerator::points(mesh::Coordinates &coordinates)
{
	_generator->points(coordinates);
}

void MeshGenerator::elements(std::vector<mesh::Element*> &elements, std::vector<eslocal> &parts)
{
	_generator->elements(elements, parts);
}

void MeshGenerator::fixPoints(std::vector<std::vector<eslocal> > &fixPoints)
{
	_generator->fixPoints(fixPoints);
}

void MeshGenerator::boundaryConditions(mesh::Coordinates &coordinates)
{
	_generator->boundaryConditions(coordinates);
}

void MeshGenerator::corners(mesh::Boundaries &boundaries)
{
	_generator->corners(boundaries);
}

void MeshGenerator::clusterBoundaries(mesh::Boundaries &boundaries)
{
	_generator->clusterBoundaries(boundaries);
}
