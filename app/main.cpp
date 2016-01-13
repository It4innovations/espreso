
#include "factory/factory.h"

int main(int argc, char **argv)
{
	MPI_Init(&argc, &argv);

	Factory factory(argc, argv);
	factory.solve(1);
	factory.store("mesh");

	MPI_Barrier(MPI_COMM_WORLD);
	MPI_Finalize();
	return 0;
}


