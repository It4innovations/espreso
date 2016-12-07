
#include "mpi.h"

#include "esinput.h"
#include "esoutput.h"
#include "esmesh.h"
#include "factory/factory.h"
#include "../config/description.h"

using namespace espreso;

int main(int argc, char** argv)
{
	if (argc < 4) {
		ESINFO(GLOBAL_ERROR) << "Specify parameters: INPUT_LOCATION  OUTPUT_LOCATION  [ NUMBER_OF_PARTS ]";
	}

	MPI_Init(&argc, &argv);

	GlobalConfiguration configuration(&argc, &argv);

	if (environment->MPIsize > 1) {
		ESINFO(GLOBAL_ERROR) << "Not implemented decomposition of ESDATA";
		configuration.input = INPUT::ESDATA;
		configuration.esdata.domains = 1;
		configuration.esdata.path = argv[1];
	} else {
		configuration.input = INPUT::WORKBENCH;
		configuration.workbench.domains = 1;
		configuration.workbench.path = argv[1];
	}

	Factory factory(configuration);
	std::cout << "Mesh loaded\n";

	for (size_t i = 2; i < argc; i++) {
		int parts = atoi(argv[i]);
		std::stringstream ss;
		ss << argv[2] << parts * environment->MPIsize;

		factory.mesh.partitiate(parts);
		std::cout << "Mesh partitiated to " << parts * environment->MPIsize << " parts\n";
		std::vector<size_t> sizes(factory.mesh.parts());
		for (size_t p = 0; p < factory.mesh.parts(); p++) {
			sizes[p] = factory.mesh.coordinates().localSize(p);
		}
		std::cout << "Nodes in subdomains: " << Info::averageValues(sizes) << "\n";
		store::Esdata::mesh(factory.mesh, ss.str());
		std::cout << "Mesh partitiated to " << parts * environment->MPIsize << " parts saved\n";
	}

	MPI_Finalize();
}


