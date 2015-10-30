#include "mpi.h"

#include "esassemblers.h"
#include "esinput.h"
#include "esoutput.h"

#include "essolver.h"
#include "esmesh.h"


int main(int argc, char** argv)
{
	int MPIrank, MPIsize;

	MPI_Init(&argc, &argv);

	MPI_Comm_rank(MPI_COMM_WORLD, &MPIrank);
	MPI_Comm_size(MPI_COMM_WORLD, &MPIsize);

	mesh::Mesh m(MPIrank, MPIsize);
	m.load(mesh::ESPRESO_INPUT, argc, argv);
	m.partitiate(16, 4);

	assembler::LinearElasticity<assembler::BEM> fem(m);


	fem.init();
	fem.solve();
	fem.post_solve_update();
	fem.finalize();

	MPI_Finalize();
}


