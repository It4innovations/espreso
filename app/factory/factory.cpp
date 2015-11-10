
#include "factory.h"

Factory::Factory(int *argc, char ***argv): _mesh(NULL), _surface(NULL)
{
	int MPIrank, MPIsize;

	MPI_Comm_rank(MPI_COMM_WORLD, &MPIrank);
	MPI_Comm_size(MPI_COMM_WORLD, &MPIsize);

	_mesh = new mesh::Mesh(MPIrank, MPIsize);
	_mesh->load(mesh::MESH_GENERATOR, *argc, *argv);

	switch (esconfig::discretization) {
		case esconfig::FEM: {
			assembler::FEM fem(*_mesh);
			_assembler = new assembler::LinearElasticity<assembler::FEM>(fem);
			break;
		}

		case esconfig::BEM: {
			_surface = new mesh::SurfaceMesh(*_mesh);
			assembler::BEM bem(*_mesh, *_surface);
			_assembler = new assembler::LinearElasticity<assembler::BEM>(bem);
			break;
		}
	}

}

Factory::~Factory()
{
	delete _assembler;
	if (_mesh != NULL) {
		delete _mesh;
	}
	if (_surface != NULL) {
		delete _surface;
	}
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



