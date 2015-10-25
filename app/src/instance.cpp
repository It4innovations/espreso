
#include "instance.h"

Instance::Instance(int argc, char** argv, int rank, int size)
	: _mesh(rank, size), _rank(rank), _size(size), _surfaceMesh(rank, size)
{

	//	_mesh.load(mesh::MESH_GENERATOR, argc, argv);

//	_mesh.load(mesh::ANSYS, argc, argv);
//	_mesh.partitiate(atoi(argv[2]), atoi(argv[3])); // for BEM


	_mesh.load(mesh::ESPRESO_INPUT, argc, argv);
	_mesh.partitiate(atoi(argv[2]), atoi(argv[3]));

	//_mesh.partitiate(2,0); // for BEM
	//_mesh.partitiate(2,8); // for FEM

}



