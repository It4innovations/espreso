#include "mpi.h"

#include "../optimization/algorithm.h"


using namespace espreso;

int main(int argc, char **argv)
{
	MPI_Init(&argc, &argv);

    int rank;
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);

    if (rank == 0)
    {
        BoothProblem bp;
        PSOAlgorithm pso(bp);
        pso.run();
    }   

    MPI_Barrier(MPI_COMM_WORLD);
	MPI_Finalize();

	return 0;
}