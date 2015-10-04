
#include "generator.h"
#include "factory/factory.h"

MeshGenerator::MeshGenerator(int argc, char** argv)
{
	_generator = MeshFactory::create(argc, argv);
}

void MeshGenerator::points(std::vector<mesh::Point> &data)
{

}

