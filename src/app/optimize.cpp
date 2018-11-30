#include "mpi.h"

#include <iostream>

#include "../optimization/optimizer.h"

using namespace espreso;

int main(int argc, char **argv)
{
	MPI_Init(&argc, &argv);

    int rank;
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);

    if (rank == 0)
    {
        Optimizer opt;
        for (int i = 0; i < 100; i++)
        { opt.set(); opt.run([&] () { std::cout << "run" << std::endl; }); }
    }   

    MPI_Barrier(MPI_COMM_WORLD);
	MPI_Finalize();

	return 0;
}