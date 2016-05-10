
#include "mpi.h"

#include "esinput.h"
#include "esoutput.h"
#include "esmesh.h"

using namespace espreso;

int main(int argc, char** argv)
{

	MPI_Init(&argc, &argv);
	Options options(&argc, &argv);

	if (options.nameless.size() < 2) {
		ESINFO(ERROR) << "Specify parameters: ANSYS_FILE OUTPUT_LOCATION [ NUMBER_OF_PARTS ]";
	}

	config::mesh::subdomains = 1;
	config::mesh::fixPoints = 0;

	// turn off compute corners
	config::solver::FETI_METHOD = config::TOTAL_FETI;

	Mesh m;
	input::AnsysWorkbench::load(m, options, 0, 1);
	std::cout << "Mesh loaded\n";

	for (size_t i = 1; i < options.nameless.size(); i++) {
		m.partitiate(atoi(options.nameless[i].c_str()));
		std::cout << "Mesh partitiated to " << options.nameless[i] << " parts\n";
		output::Esdata data(m, (options.nameless[0] + options.nameless[i]).c_str());
		data.store(1, 1);
		std::cout << "Mesh partitiated to " << options.nameless[i] << " parts saved\n";
	}

	MPI_Finalize();
}


