#include "mpi.h"

#include <iostream>
#include <fstream>
#include <sstream>

#include "../config/export/json.h"

using namespace espreso;

int main(int argc, char **argv)
{
	MPI_Init(&argc, &argv);

    int rank;
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);

    ECF ecf;
    if (rank == 0)
    {
        ECFJSONExport json(&ecf, std::cout);
        json.exportToStream();

		// ecf.forEachParameters([&] (const ECFParameter *parameter) {
		// 	std::cout << parameter->name << std::endl;
		// 	if (parameter->metadata.condition->isset()) {
		// 		std::cout << parameter->name << ": ";
		// 		ecf.forEachParameters([&] (const ECFParameter *conditionalParameter) {
		// 			if (parameter->metadata.condition->match(conditionalParameter->data())) {
		// 				std::cout << parameter->metadata.condition->compose(conditionalParameter);
		// 			}
		// 		}, false, true);
		// 		std::cout << "\n";
		// 	}
		// }, false, true);
    }   

    MPI_Barrier(MPI_COMM_WORLD);
	MPI_Finalize();

	return 0;
}



