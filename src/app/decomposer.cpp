
#include "mpi.h"

#include "esinput.h"
#include "esoutput.h"
#include "esmesh.h"
#include "factory/factory.h"

using namespace espreso;

int main(int argc, char** argv)
{

	MPI_Init(&argc, &argv);

	if (config::env::MPIsize > 1) {
		config::mesh::INPUT = config::mesh::ESDATA;
	} else {
		config::mesh::INPUT = config::mesh::WORKBENCH;
	}

	Configuration configuration = ParametersReader::arguments(&argc, &argv);

	if (configuration.nameless.size() < 2) {
		ESINFO(ERROR) << "Specify parameters: INPUT_LOCATION  OUTPUT_LOCATION  [ NUMBER_OF_PARTS ]";
	}

	config::mesh::SUBDOMAINS = 1;
	config::mesh::FIX_POINTS = 0;
	// turn off compute corners
	config::solver::FETI_METHOD = config::TOTAL_FETI;

	Factory factory(configuration);
	std::cout << "Mesh loaded\n";

	for (size_t i = 1; i < configuration.nameless.size(); i++) {
		int parts = atoi(configuration.nameless[i].c_str());
		std::stringstream ss;
		ss << configuration.nameless[0] << parts * config::env::MPIsize;

		factory.mesh()->partitiate(parts);
		std::cout << "Mesh partitiated to " << parts * config::env::MPIsize << " parts\n";
		output::Esdata data(*factory.mesh(), ss.str());
		std::vector<size_t> sizes(factory.mesh()->parts());
		for (size_t p = 0; p < factory.mesh()->parts(); p++) {
			sizes[p] = factory.mesh()->getPartNodesCount(p);
		}
		std::cout << Info::averageValues(sizes) << "\n";
		data.store(1, 1);
		std::cout << "Mesh partitiated to " << parts * config::env::MPIsize << " parts saved\n";
	}

	MPI_Finalize();
}


