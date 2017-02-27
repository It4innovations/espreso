
#include "mpi.h"

#include "../configuration/globalconfiguration.h"
#include "../mesh/structures/mesh.h"
#include "../mesh/structures/coordinates.h"
#include "factory/factory.h"
#include "../output/espreso/espresobinaryformat.h"
#include "../basis/logging/logging.hpp"

using namespace espreso;

int main(int argc, char** argv)
{
	MPI_Init(&argc, &argv);
	GlobalConfiguration configuration(&argc, &argv);

	Factory factory(configuration);
	std::stringstream decomposition(configuration.decomposer.parts);

	size_t parts;
	while (decomposition >> parts) {
		std::stringstream path;
		path << configuration.decomposer.prefix << parts * environment->MPIsize;

		factory.mesh->partitiate(parts);
		ESINFO(ALWAYS) << "Mesh partitiated to " << parts * environment->MPIsize << " parts";
		std::vector<size_t> sizes(factory.mesh->parts());
		for (size_t p = 0; p < factory.mesh->parts(); p++) {
			sizes[p] = factory.mesh->coordinates().localSize(p);
		}
		ESINFO(ALWAYS) << "Nodes in domains: " << Info::averageValues(sizes);
		store::ESPRESOBinaryFormat::store(*factory.mesh, path.str());
		ESINFO(ALWAYS) << "Mesh partitiated to " << parts * environment->MPIsize << " parts saved";
	}

	factory.finalize();
	MPI_Finalize();
}


