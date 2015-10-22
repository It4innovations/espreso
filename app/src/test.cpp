#include "mpi.h"

#include "esinput.h"
#include "esoutput.h"

#include "essolver.h"
#include "esmesh.h"
#include "esphysics.h"
#include "escomposer.h"


int main(int argc, char** argv)
{
	int MPIrank, MPIsize;

	MPI_Init(&argc, &argv);

	MPI_Comm_rank(MPI_COMM_WORLD, &MPIrank);
	MPI_Comm_size(MPI_COMM_WORLD, &MPIsize);

	mesh::Mesh m(MPIrank, MPIsize);
	m.load(mesh::MESH_GENERATOR, argc, argv);

	composer::FEM<physics::LinearElasticity> fem(m);

	SparseMatrix sm;
	fem.element(sm, 0);

	MPI_Finalize();
}


