
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

static void buildSequential(OrderedMeshDatabase &database, Mesh &mesh)
{
	OrderedMesh _mesh;
}

void build(OrderedMeshDatabase &database, Mesh &mesh)
{
	if (info::mpi::size == 1) {
		buildSequential(database, mesh);
		return;
	}

	eslog::startln("BUILDER: BUILD SCATTERED MESH", "BUILDER");

	OrderedMesh ordered;

	balance(database, ordered);
	eslog::checkpointln("BUILDER: DATA BALANCED");

	Communication::allReduce(&ordered.dimension, NULL, 1, MPI_INT, MPI_MAX); // we reduce dimension here in order to be able measure inbalances in the 'balance' function
	HilbertCurve<esfloat> sfc(ordered.dimension, SFCDEPTH, ordered.coordinates.size(), ordered.coordinates.data());

	eslog::checkpointln("BUILDER: SFC BUILT");

	ClusteredMesh clustered;

	assignBuckets(ordered, sfc, clustered);
	eslog::checkpointln("BUILDER: SFC BUCKETS ASSIGNED");

	clusterize(ordered, clustered);
	eslog::checkpointln("BUILDER: MESH CLUSTERIZED");

	computeSFCNeighbors(sfc, clustered);
	eslog::checkpointln("BUILDER: NEIGHBORS APPROXIMATED");

	mergeDuplicatedNodes(sfc, clustered);
	eslog::checkpointln("BUILDER: DUPLICATED NODES FOUND");

	groupElementTypes(clustered);
	eslog::checkpointln("BUILDER: ELEMENTS GROUPED");

	linkup(clustered);
	eslog::checkpointln("BUILDER: LINKED UP");

	MPI_Finalize();
	exit(0);

////	balance();
////	eslog::checkpointln("BUILDER: DATA BALANCED");
//
//	assignRegions(_meshData.eregions, _meshData.eIDs, _eDistribution, _eregsize, _eregions);
//	assignRegions(_meshData.nregions, _meshData.nIDs, _nDistribution, _nregsize, _nregions);
//	eslog::checkpointln("BUILDER: REGION ASSIGNED");
//
////	reindexRegions();
//
//	assignNBuckets();
//	eslog::checkpointln("BUILDER: NODES BUCKETS COMPUTED");
//
//	assignEBuckets();
//	eslog::checkpointln("BUILDER: ELEMENTS BUCKETS ASSIGNED");
//
//	clusterize();
//	eslog::checkpointln("BUILDER: ELEMENTS CLUSTERED");
//
//	computeSFCNeighbors();
//	eslog::checkpointln("BUILDER: NEIGHBORS APPROXIMATED");
//
//	if (_meshData.removeDuplicates) {
//		mergeDuplicatedNodes();
//		eslog::checkpointln("BUILDER: DUPLICATED NODES MERGED");
//	}
//
//	sortElements();
//	eslog::checkpointln("BUILDER: ELEMENTS SORTED");
//
//	linkup();
//	fillNeighbors();
//	eslog::checkpointln("BUILDER: LINKED UP");
//
//	sortNodes();
//	eslog::checkpointln("BUILDER: NODES SORTED");
//
//	reindexElementNodes();
//	eslog::checkpointln("BUILDER: ELEMENTS NODES REINDEXED");
//
//	if (_meshData.removeDuplicates) {
//		removeDuplicateElements();
//		eslog::checkpointln("BUILDER: DUPLICATED ELEMENTS REMOVED");
//	}
//
//	fillNodes();
//	eslog::checkpointln("BUILDER: NODES FILLED");
//
//	fillElements();
//	eslog::checkpointln("BUILDER: ELEMENTS SORTED");
//
//	if (info::mesh->nodes->elements == NULL) {
//		mesh::linkNodesAndElements(info::mesh->elements, info::mesh->nodes, info::mesh->neighbors);
//	}
//
//	exchangeBoundary();
//	eslog::checkpointln("BUILDER: BOUNDARY EXCHANGED");
//
//	fillRegions(_meshData.eregions, _eregsize, _eregions);
//	fillRegions(_meshData.nregions, _nregsize, _nregions);
//	fillElementRegions();
//	fillBoundaryRegions();
//	fillNodeRegions();
//	eslog::checkpointln("BUILDER: REGIONS FILLED");
//
//	reindexBoundaryNodes();
//	eslog::endln("BUILDER: BOUNDARY NODES REINDEXED");
}

}
}
