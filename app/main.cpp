
#include "factory/factory.h"

int main(int argc, char **argv)
{
	MPI_Init(&argc, &argv);

	//Factory app(&argc, &argv);
	//app.solve(1);


	int MPIrank, MPIsize;

	MPI_Comm_rank(MPI_COMM_WORLD, &MPIrank);
	MPI_Comm_size(MPI_COMM_WORLD, &MPIsize);

	mesh::Mesh mesh(MPIrank, MPIsize);
	mesh.load(mesh::MESH_GENERATOR, argc, argv);

	assembler::FEM fem(mesh);
	assembler::LinearElasticity<assembler::FEM> app(fem);

	app.init();
	app.solve();
	app.post_solve_update();
	app.finalize();

	MPI_Barrier(MPI_COMM_WORLD);
	MPI_Finalize();
	return 0;
}


