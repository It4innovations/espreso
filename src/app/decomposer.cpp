
#include "mpi.h"

#include "../mesh/structures/mesh.h"
#include "../mesh/structures/coordinates.h"
#include "../output/esdata/esdata.h"
#include "factory/factory.h"
#include "../config/globalconfiguration.h"

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
		store::Esdata::mesh(*factory.mesh, path.str());
		ESINFO(ALWAYS) << "Mesh partitiated to " << parts * environment->MPIsize << " parts saved";
	}

	MPI_Finalize();
}


