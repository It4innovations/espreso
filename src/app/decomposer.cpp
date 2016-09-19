
#include "mpi.h"

#include "esinput.h"
#include "esoutput.h"
#include "esmesh.h"
#include "factory/factory.h"

using namespace espreso;

int main(int argc, char** argv)
{

	MPI_Init(&argc, &argv);
	Options options(&argc, &argv);

	if (options.nameless.size() < 2) {
		ESINFO(ERROR) << "Specify parameters: INPUT_LOCATION  OUTPUT_LOCATION  [ NUMBER_OF_PARTS ]";
	}

	if (config::env::MPIsize > 1) {
		ESINFO(GLOBAL_ERROR) << "Not implemented decomposition of ESDATA";
		config::mesh::INPUT = config::mesh::INPUTalternative::ESDATA;
	} else {
		config::mesh::INPUT = config::mesh::INPUTalternative::WORKBENCH;
	}
	config::mesh::SUBDOMAINS = 1;
	config::info::VERBOSE_LEVEL = 2;
	config::info::MEASURE_LEVEL = 2;

	Factory factory(configuration);
	std::cout << "Mesh loaded\n";

	for (size_t i = 1; i < options.nameless.size(); i++) {
		int parts = atoi(options.nameless[i].c_str());
		std::stringstream ss;
		ss << options.nameless[0] << parts * config::env::MPIsize;

		factory.mesh.partitiate(parts);
		std::cout << "Mesh partitiated to " << parts * config::env::MPIsize << " parts\n";
		std::vector<size_t> sizes(factory.mesh.parts());
		for (size_t p = 0; p < factory.mesh.parts(); p++) {
			sizes[p] = factory.mesh.coordinates().localSize(p);
		}
		std::cout << "Nodes in subdomains: " << Info::averageValues(sizes) << "\n";
		output::Esdata::mesh(factory.mesh, ss.str());
		std::cout << "Mesh partitiated to " << parts * config::env::MPIsize << " parts saved\n";
	}

	MPI_Finalize();
}


