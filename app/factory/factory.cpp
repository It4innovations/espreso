
#include "factory.h"

Factory::Factory(int argc, char **argv)
{
	int MPIrank, MPIsize;

	MPI_Comm_rank(MPI_COMM_WORLD, &MPIrank);
	MPI_Comm_size(MPI_COMM_WORLD, &MPIsize);

	mesh::Mesh mesh(MPIrank, MPIsize);
	mesh.load(mesh::MESH_GENERATOR, argc, argv);

	assembler::FEM fem(mesh);
	_assembler = new assembler::LinearElasticity<assembler::FEM>(fem);
}

Factory::~Factory()
{
	delete _assembler;
}

void Factory::solve(eslocal steps)
{
	_assembler->init();

	for (int i = 0; i < steps; i++) {
		_assembler->pre_solve_update();
		_assembler->solve();
		_assembler->post_solve_update();
	}

	_assembler->finalize();
}



