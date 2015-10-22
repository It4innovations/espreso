
#include "instance.h"

Instance::Instance(int argc, char** argv, int rank, int size)
	: _mesh(rank, size), _rank(rank), _size(size), _surfaceMesh(rank, size)
{
	_mesh.load(mesh::MESH_GENERATOR, argc, argv);
	//_mesh.store(mesh::VTK, "sphere", 0.9, 0.9);

//	_mesh.load(mesh::ESPRESO_INPUT, argc, argv);
//	_mesh.partitiate(2,0); // for BEM
//	//_mesh.partitiate(32,8); // for FEM

}



