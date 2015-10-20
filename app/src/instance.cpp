
#include "instance.h"

Instance::Instance(int argc, char** argv, int rank, int size)
	: _mesh(rank, size), _rank(rank), _size(size)
{
	_mesh.load(mesh::MESH_GENERATOR, argc, argv);
}



