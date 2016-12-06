
#include "mpi.h"

#include "esinput.h"
#include "esoutput.h"
#include "esmesh.h"
#include "factory/factory.h"

using namespace espreso;

int main(int argc, char** argv)
{
	MPI_Init(&argc, &argv);

	ArgsConfiguration configuration = ParametersReader::fromArguments(&argc, &argv);

	config::mesh::INPUT = config::mesh::INPUTalternative::GENERATOR;
	config::mesh::SUBDOMAINS = 1;
	config::mesh::FIX_POINTS = 0;
	config::solver::FETI_METHOD = config::solver::FETI_METHODalternative::TOTAL_FETI;
	config::info::VERBOSE_LEVEL = 3;
	config::info::MEASURE_LEVEL = 2;

	Mesh mesh;
	input::CubeSettings cube(configuration, 0, 1);
	input::CubeGenerator<input::Hexahedron8>::load(mesh, cube);
	config::info::VERBOSE_LEVEL = 2;
	config::info::MEASURE_LEVEL = 2;

	std::cout << "Mesh loaded\n";

	std::cout << "Mesh partitiated to " << mesh.parts() * environment->MPIsize << " parts\n";
	std::vector<size_t> sizes(mesh.parts());
	for (size_t p = 0; p < mesh.parts(); p++) {
		sizes[p] = mesh.coordinates().localSize(p);
	}
	std::cout << "Nodes in subdomains: " << Info::averageValues(sizes) << "\n";
	std::stringstream ss;
	ss << configuration.nameless.back() << mesh.parts();
	output::Esdata::mesh(mesh, ss.str());
	std::cout << "Mesh saved\n";

	MPI_Finalize();
}






