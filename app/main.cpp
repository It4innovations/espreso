
#include "factory/factory.h"

int main(int argc, char **argv)
{
	MPI_Init(&argc, &argv);

	Factory app(argc, argv);

	app.solve(1);

	MPI_Barrier(MPI_COMM_WORLD);
	MPI_Finalize();
	return 0;
}


