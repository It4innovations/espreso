
#include "factory/factory.h"

int main(int argc, char **argv)
{
	MPI_Init(&argc, &argv);

	Options options(&argc, &argv);

	Factory factory(options);
	factory.solve(1);
	factory.store("mesh");

	MPI_Barrier(MPI_COMM_WORLD);
	ESLOG(eslog::DURATION) << "ESPRESO overall time";
	MPI_Finalize();
	return 0;
}


