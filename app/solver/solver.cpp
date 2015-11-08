
#include "solver.h"

Espreso::Espreso(int argc, char **argv)
{
	int MPIrank, MPIsize;

	MPI_Init(&argc, &argv);

	MPI_Comm_rank(MPI_COMM_WORLD, &MPIrank);
	MPI_Comm_size(MPI_COMM_WORLD, &MPIsize);

	mesh::Mesh mesh(MPIrank, MPIsize);
	mesh.load(mesh::MESH_GENERATOR, argc, argv);

	assembler::FEM fem(mesh);
	assembler::LinearElasticity<assembler::FEM> solver(fem);

	solver.init();
	solver.solve();
	solver.post_solve_update();
	solver.finalize();

	MPI_Finalize();
}



