
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

void build(InputMesh<OrderedUniqueNodes, OrderedUniqueFaces, OrderedRegions> &input, Mesh &mesh)
{
	eslog::startln("BUILDER: PROCESS FACED MESH", "BUILDER");
	mesh.dimension = 3;
	if (info::mpi::size == 1) {
		TemporalSequentialMesh<MergedNodes, OrderedFacesBalanced> grouped;
		TemporalSequentialMesh<MergedNodes, MergedElements> prepared;

		initialize(input, grouped);
		eslog::checkpointln("BUILDER: DATA INITIALIZED");

		buildElementsFromFaces(grouped, prepared);
		eslog::checkpointln("BUILDER: ELEMENTS CREATED");

		fillMesh(prepared, *input.regions, mesh);
	} else {
		TemporalMesh<OrderedNodesBalanced, OrderedFacesBalanced> grouped;
		TemporalMesh<OrderedNodesBalanced, OrderedElementsBalanced> ordered;
		TemporalMesh<MergedNodes, MergedElements> clustered;
		TemporalMesh<LinkedNodes, MergedElements> linked;

		balance(input, grouped);
		eslog::checkpointln("BUILDER: DATA BALANCED");

		buildElementsFromFaces(grouped, ordered);
		eslog::checkpointln("BUILDER: ELEMENTS CREATED");

		ivector<esint> nbuckets, ebuckets, splitters;
		HilbertCurve<esfloat> sfc(mesh.dimension, SFCDEPTH, ordered.nodes->coordinates.size(), ordered.nodes->coordinates.data());
		eslog::checkpointln("BUILDER: SFC BUILT");

		assignBuckets(ordered, sfc, nbuckets, ebuckets);
		eslog::checkpointln("BUILDER: SFC BUCKETS ASSIGNED");

		clusterize(ordered, nbuckets, ebuckets, sfc.buckets(sfc.depth), clustered, splitters);
		utils::clearVectors(nbuckets, ebuckets);
		eslog::checkpointln("BUILDER: MESH CLUSTERIZED");

		computeSFCNeighbors(sfc, splitters, linked.nodes->neighbors); // neighbors are approximated here
		eslog::checkpointln("BUILDER: NEIGHBORS APPROXIMATED");

		linkup(clustered, linked);
		reindexToLocal(linked);
		eslog::checkpointln("BUILDER: LINKED UP");

		fillMesh(linked, *input.regions, mesh);
	}
	eslog::endln("BUILDER: MESH BUILT");
}

}
}
