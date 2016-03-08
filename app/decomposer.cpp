
#include "mpi.h"

#include "esinput.h"
#include "esoutput.h"
#include "esmesh.h"

int main(int argc, char** argv)
{
	if (argc < 4) {
		ESLOG(eslog::ERROR) << "Specify parameters: ANSYS_FILE OUTPUT_LOCATION [ NUMBER_OF_PARTS ]";
	}


	MPI_Init(&argc, &argv);

	Options options(&argc, &argv);

	esconfig::mesh::subdomains = 1;
	esconfig::mesh::fixPoints = 0;

	// turn off compute corners
	esconfig::solver::FETI_METHOD = esconfig::TOTAL_FETI;

	mesh::Mesh m;
	esinput::AnsysWorkbench loader(options, 0, 1);
	loader.load(m);
	std::cout << "Mesh loaded\n";

	for (size_t i = 1; i < options.nameless.size(); i++) {
		m.partitiate(atoi(options.nameless[i].c_str()));
		std::cout << "Mesh partitiated to " << options.nameless[i] << " parts\n";
		esoutput::Esdata data(m, (options.nameless[0] + options.nameless[i]).c_str());
		data.store(1, 1);
		std::cout << "Mesh partitiated to " << options.nameless[i] << " parts saved\n";
	}

	MPI_Finalize();
}


