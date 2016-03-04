
#include "mpi.h"

#include "esinput.h"
#include "esoutput.h"
#include "esmesh.h"

int main(int argc, char** argv)
{
	if (argc < 4) {
		ESLOG(eslog::ERROR) << "Specify parameters: ANSYS_FILE NUMBER_OF_PARTS OUTPUT_LOCATION";
	}

	MPI_Init(&argc, &argv);

	Options options(&argc, &argv);

	esconfig::mesh::subdomains = atoi(options.nameless[0].c_str());
	esconfig::mesh::fixPoints = 0;

	// turn off compute corners
	esconfig::solver::FETI_METHOD = esconfig::TOTAL_FETI;

	mesh::Mesh m;
	esinput::AnsysWorkbench loader(options, 0, 1);
	loader.load(m);

	esoutput::VTK_Full vtk(m, "decomposer");
	vtk.store(0.9, 0.9);

	esoutput::Esdata data(m, options.nameless[1].c_str());
	data.store(1, 1);

	MPI_Finalize();
}


