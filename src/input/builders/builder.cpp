
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

namespace espreso {
namespace builder {

//static void buildSequential(OrderedMeshDatabase &database, Mesh &mesh)
//{
//	OrderedMesh _mesh;
//}

void build(InputMesh<OrderedNodes, OrderedElements, OrderedRegions> &input, Mesh &mesh)
{
	if (info::mpi::size == 1) {
//		buildSequential(database, mesh);
		return;
	}

	eslog::startln("BUILDER: PROCESS ORDERED MESH", "BUILDER");

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
	ordered.clear();
	utils::clearVectors(nbuckets, ebuckets);
	eslog::checkpointln("BUILDER: MESH CLUSTERIZED");

	computeSFCNeighbors(sfc, clustered, splitters, linked.nodes->neighbors); // neighbors are approximated here
	eslog::checkpointln("BUILDER: NEIGHBORS APPROXIMATED");

	searchDuplicatedNodes(sfc, splitters, linked.nodes->neighbors, clustered.nodes, merged.nodes);
	delete merged.elements; merged.elements = clustered.elements; clustered.elements = nullptr;
	clustered.clear();
	utils::clearVector(splitters);
	eslog::checkpointln("BUILDER: DUPLICATED NODES FOUND");

	linkup(merged, linked);
	merged.clear();
	eslog::checkpointln("BUILDER: LINKED UP");

	searchParentAndDuplicatedElements(linked, prepared, mesh.dimension);
	linked.clear();
	eslog::checkpointln("BUILDER: DUPLICATED ELEMENTS FOUND");

	fillMesh(prepared, *input.regions, mesh);
	prepared.clear();
	eslog::endln("BUILDER: MESH BUILT");
}

}
}
