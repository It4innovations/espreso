
#include "builder.h"
#include "builder.utils.h"

#include "basis/containers/serializededata.h"
#include "basis/sfc/hilbertcurve.h"
#include "basis/utilities/packing.h"
#include "basis/utilities/utils.h"
#include "esinfo/eslog.h"
#include "esinfo/envinfo.h"
#include "wrappers/mpi/communication.h"

#include <numeric>
#include <algorithm>

#include "basis/utilities/print.h"

namespace espreso {
namespace builder {

void build(InputMesh<OrderedNodes, OrderedElements, OrderedRegions> &input, Mesh &mesh)
{
	eslog::startln("BUILDER: PROCESS ORDERED MESH", "BUILDER");
	if (info::mpi::size == 1) {
		TemporalSequentialMesh<ClusteredNodes, ClusteredElements> clustered;
		TemporalSequentialMesh<MergedNodes, ClusteredElements> merged;
		TemporalSequentialMesh<MergedNodes, MergedElements> prepared;

		initialize(input, clustered, mesh.dimension);
		eslog::checkpointln("BUILDER: DATA INITIALIZED");

		searchDuplicatedNodes(clustered, merged);
		eslog::checkpointln("BUILDER: DUPLICATED NODES FOUND");

		searchDuplicatedElements(merged, prepared, mesh.dimension);
		eslog::checkpointln("BUILDER: DUPLICATED ELEMENTS FOUND");

		fillMesh(prepared, *input.regions, mesh);
	} else {
		TemporalMesh<OrderedNodesBalanced, OrderedElementsBalanced> ordered;
		TemporalMesh<ClusteredNodes, ClusteredElements> clustered;
		TemporalMesh<MergedNodes, ClusteredElements> merged;
		TemporalMesh<LinkedNodes, ClusteredElements> linked;
		TemporalMesh<LinkedNodes, MergedElements> prepared;

		ivector<esint> nbuckets, ebuckets, splitters;

		balance(input, ordered, mesh.dimension);
		eslog::checkpointln("BUILDER: DATA BALANCED");

		Communication::allReduce(&mesh.dimension, NULL, 1, MPI_INT, MPI_MAX); // we reduce dimension here in order to be able measure inbalances in the 'balance' function
		HilbertCurve<esfloat> sfc(mesh.dimension, SFCDEPTH, ordered.nodes->coordinates.size(), ordered.nodes->coordinates.data());
		eslog::checkpointln("BUILDER: SFC BUILT");

		assignBuckets(ordered, sfc, nbuckets, ebuckets);
		eslog::checkpointln("BUILDER: SFC BUCKETS ASSIGNED");

		clusterize(ordered, nbuckets, ebuckets, sfc.buckets(sfc.depth), clustered, splitters);
		utils::clearVectors(nbuckets, ebuckets);
		eslog::checkpointln("BUILDER: MESH CLUSTERIZED");

		computeSFCNeighbors(sfc, splitters, linked.nodes->neighbors); // neighbors are approximated here
		eslog::checkpointln("BUILDER: NEIGHBORS APPROXIMATED");

		searchDuplicatedNodes(sfc, splitters, linked.nodes->neighbors, clustered, merged);
		utils::clearVector(splitters);
		eslog::checkpointln("BUILDER: DUPLICATED NODES FOUND");

		linkup(merged, linked);
		eslog::checkpointln("BUILDER: LINKED UP");

		searchParentAndDuplicatedElements(linked, prepared, mesh.dimension);
		eslog::checkpointln("BUILDER: DUPLICATED ELEMENTS FOUND");

		fillMesh(prepared, *input.regions, mesh);
	}
	eslog::endln("BUILDER: MESH BUILT");
}

void build(InputMesh<OrderedUniqueNodes, OrderedUniqueFaces, OrderedUniqueFacesRegions> &input, Mesh &mesh)
{
	eslog::startln("BUILDER: PROCESS FACED MESH", "BUILDER");
	if (info::mpi::size == 1) {
		TemporalSequentialMesh<ClusteredNodes, ClusteredElements> prepared;

		initialize(input, prepared, mesh.dimension);
		eslog::checkpointln("BUILDER: DATA INITIALIZED");

		input.regions->files.wait();
//		fillMesh(prepared, *input.regions, mesh);
	} else {

	}
	eslog::endln("BUILDER: MESH BUILT");

	MPI_Finalize();
	exit(0);
}

}
}
