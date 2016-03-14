
#include "factory/factory.h"

using namespace espreso;

int main(int argc, char **argv)
{
	MPI_Init(&argc, &argv);

	Options options(&argc, &argv);

	Factory factory(options);
	factory.solve(1);
	factory.store("result");

	MPI_Barrier(MPI_COMM_WORLD);
	MPI_Finalize();
	return 0;
}


