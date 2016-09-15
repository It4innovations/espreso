
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

//	if (config::env::MPIsize > 1) {
//		config::mesh::INPUT = config::mesh::INPUTalternative::ESDATA;
//	} else {
//		config::mesh::INPUT = config::mesh::INPUTalternative::WORKBENCH;
//	}
	config::mesh::INPUT = config::mesh::INPUTalternative::GENERATOR;
	config::mesh::SUBDOMAINS = 1;
	config::mesh::FIX_POINTS = 0;
	config::solver::FETI_METHOD = config::solver::FETI_METHODalternative::TOTAL_FETI;
	config::info::VERBOSE_LEVEL = 2;
	config::info::MEASURE_LEVEL = 2;

	//Factory factory(configuration);
	Mesh mesh;

	input::CubeSettings cube(0, 1);
	cube.nodes["BOTTOM"] = Interval(0, cube.problemLength[0], 0, cube.problemLength[1], 0, 0);
	cube.properties["DIRICHLET"]["BOTTOM"] = "x: 0, y: 0, z: 0";
	cube.subdomainsInCluster[0] = 2;
	cube.subdomainsInCluster[1] = 2;
	cube.subdomainsInCluster[2] = 2;
	cube.elementsInSubdomain[0] = 8;
	cube.elementsInSubdomain[1] = 8;
	cube.elementsInSubdomain[2] = 8;
	input::CubeGenerator<input::Hexahedron8>::load(mesh, cube);

	std::cout << "Mesh loaded\n";

	for (size_t i = 1; i < options.nameless.size(); i++) {
		int parts = atoi(options.nameless[i].c_str());
		std::stringstream ss;
		ss << options.nameless[0] << parts * config::env::MPIsize;

		mesh.partitiate(parts);
		std::cout << "Mesh partitiated to " << parts * config::env::MPIsize << " parts\n";
		std::vector<size_t> sizes(mesh.parts());
		for (size_t p = 0; p < mesh.parts(); p++) {
			sizes[p] = mesh.coordinates().localSize(p);
		}
		std::cout << "Nodes in subdomains: " << Info::averageValues(sizes) << "\n";
		output::Esdata::mesh(mesh, ss.str());
		std::cout << "Mesh partitiated to " << parts * config::env::MPIsize << " parts saved\n";
	}

	MPI_Finalize();
}


