
#include "instance.h"

Instance::Instance(int argc, char** argv, int rank, int size)
	: _mesh(rank, size), _rank(rank), _size(size), _surfaceMesh(rank, size)
{
//	_mesh.load(mesh::MESH_GENERATOR, argc, argv);

	_mesh.load(mesh::ESPRESO_INPUT, argc, argv);
	_mesh.partitiate(32,4);

}



