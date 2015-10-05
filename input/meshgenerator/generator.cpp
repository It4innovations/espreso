
#include "generator.h"
#include "factory/factory.h"

using namespace esinput;


MeshGenerator::MeshGenerator(int argc, char** argv, int rank, int size)
{
	_generator = MeshFactory::create(argc, argv, rank, size);
}

void MeshGenerator::points(mesh::Coordinates &data)
{
	_generator->points(data);
}

void MeshGenerator::elements(std::vector<mesh::Element*> &data)
{
	_generator->elements(data);
}

