
#include "mpi.h"

#include "../configuration/globalconfiguration.h"
#include "../mesh/structures/mesh.h"
#include "../mesh/structures/coordinates.h"
#include "../input/loader.h"
#include "../basis/logging/logging.hpp"
#include "../output/datastore/espresobinaryformat.h"

using namespace espreso;

int main(int argc, char** argv)
{
	MPI_Init(&argc, &argv);
	GlobalConfiguration configuration(&argc, &argv);

	size_t parts;
	std::stringstream directoryTree(configuration.decomposer.parts);
	while (directoryTree >> parts) {
		std::stringstream path;
		path << configuration.decomposer.prefix << parts * environment->MPIsize;
		ESPRESOBinaryFormat::prepareDirectories(path.str(), parts);
	}

	Mesh mesh;
	input::Loader::load(configuration, mesh, configuration.env.MPIrank, configuration.env.MPIsize);
	std::stringstream decomposition(configuration.decomposer.parts);
	while (decomposition >> parts) {
		std::stringstream path;
		path << configuration.decomposer.prefix << parts * environment->MPIsize;

		mesh.partitiate(parts);
		ESINFO(ALWAYS) << "Mesh partitiated to " << parts * environment->MPIsize << " parts";
		std::vector<size_t> sizes(mesh.parts());
		for (size_t p = 0; p < mesh.parts(); p++) {
			sizes[p] = mesh.coordinates().localSize(p);
		}
		ESINFO(ALWAYS) << "Nodes in domains: " << Info::averageValues(sizes);
		ESPRESOBinaryFormat::store(mesh, path.str());
		ESINFO(ALWAYS) << "Mesh partitiated to " << parts * environment->MPIsize << " parts saved";
	}

	MPI_Finalize();
}


