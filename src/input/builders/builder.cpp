
#include "builder.h"
#include "builder.utils.h"

#include "basis/containers/serializededata.h"
#include "basis/sfc/hilbertcurve.h"
#include "basis/utilities/packing.h"
#include "basis/utilities/utils.h"
#include "esinfo/eslog.hpp"
#include "esinfo/envinfo.h"
#include "mesh/store/elementsregionstore.h"
#include "mesh/store/boundaryregionstore.h"
#include "wrappers/mpi/communication.h"

#include <numeric>
#include <algorithm>

namespace espreso {
namespace builder {

void buildOrderedFEM(InputMesh<OrderedNodes, OrderedElements, OrderedRegions> &input, Mesh &mesh)
{
	eslog::startln("BUILDER: PROCESS ORDERED MESH", "BUILDER");
	eslog::param("OrderedNodes", size(*input.nodes));
	eslog::param("OrderedElements", size(*input.elements));
	eslog::param("OrderedTotal", size(*input.nodes) + size(*input.elements));

	eslog::info(" ==================================== ORDERED FEM BUILDER ===================== %12.3f s\n", eslog::duration());
	eslog::info(" ============================================================================================= \n");

	if (info::mpi::size == 1) {
		TemporalSequentialMesh<ClusteredNodes, ClusteredElements> clustered;
		TemporalSequentialMesh<MergedNodes, ClusteredElements> merged;
		TemporalSequentialMesh<MergedNodes, MergedElements> prepared;

		initializeSequentialFEM(input, clustered, mesh.dimension);
		eslog::info(" == MESH DIMENSION %72d == \n", mesh.dimension);
		eslog::checkpointln("BUILDER: DATA INITIALIZED");

		searchDuplicatedNodes(clustered, merged);
		eslog::checkpointln("BUILDER: DUPLICATED NODES FOUND");

		searchDuplicatedElements(merged, prepared, mesh.dimension);
		eslog::checkpointln("BUILDER: DUPLICATED ELEMENTS FOUND");

		fillSequentialMesh(prepared, *input.regions, mesh);
	} else {
		TemporalMesh<OrderedNodesBalanced, OrderedElementsBalanced> ordered;
		TemporalMesh<ClusteredNodes, ClusteredElements> clustered;
		TemporalMesh<MergedNodes, ClusteredElements> merged;
		TemporalMesh<LinkedNodes, ClusteredElements> linked;
		TemporalMesh<LinkedNodes, MergedElements> prepared;

		ivector<esint> nbuckets, ebuckets, splitters;

		balanceFEM(input, ordered, mesh.dimension);
		eslog::checkpointln("BUILDER: DATA BALANCED");
		eslog::param("OrderedNodesBalanced", size(*ordered.nodes));
		eslog::param("OrderedElementsBalanced", size(*ordered.elements));
		eslog::param("OrderedTotalBalanced", size(*ordered.nodes) + size(*ordered.elements));

		Communication::allReduce(&mesh.dimension, NULL, 1, MPI_INT, MPI_MAX); // we reduce dimension here in order to be able measure inbalances in the 'balance' function
		eslog::info(" == MESH DIMENSION %72d == \n", mesh.dimension);
		HilbertCurve<esfloat> sfc(mesh.dimension, SFCDEPTH, ordered.nodes->coordinates.size(), ordered.nodes->coordinates.data());
		eslog::checkpointln("BUILDER: SFC BUILT");

		assignBuckets(ordered, sfc, nbuckets, ebuckets);
		eslog::checkpointln("BUILDER: SFC BUCKETS ASSIGNED");

		clusterize(ordered, nbuckets, ebuckets, sfc.buckets(sfc.depth), clustered, splitters);
		utils::clearVectors(nbuckets, ebuckets);
		eslog::checkpointln("BUILDER: MESH CLUSTERIZED");

		computeSFCNeighbors(sfc, splitters, linked.nodes->neighbors); // neighbors are approximated here
		eslog::checkpointln("BUILDER: NEIGHBORS APPROXIMATED");

		searchDuplicatedNodesWithSFC(sfc, splitters, linked.nodes->neighbors, clustered, merged);
		utils::clearVector(splitters);
		eslog::checkpointln("BUILDER: DUPLICATED NODES FOUND");

		linkup(merged, linked);
		eslog::checkpointln("BUILDER: LINKED UP");

		searchParentAndDuplicatedElements(linked, prepared, mesh.dimension);
		eslog::checkpointln("BUILDER: DUPLICATED ELEMENTS FOUND");

		fillMesh(prepared, *input.regions, mesh);
	}
	int bregions = 0, nregions = 0;
	for (size_t r = 1; r < mesh.boundaryRegions.size(); ++r) {
		if (mesh.boundaryRegions[r]->dimension) {
			++bregions;
		} else {
			++nregions;
		}
	}
	eslog::info(" ==    -   -   -   -   -   -   -   -   -   -   -   -   -   -   -   -   -   -   -   -   -    == \n");
	eslog::info(" == ELEMENTS REGIONS %70d == \n", mesh.elementsRegions.size() - 1);
	eslog::info(" == BOUNDARY REGIONS %70d == \n", bregions);
	eslog::info(" == NODES REGIONS %73d == \n", nregions);
	eslog::info(" ============================================================================================= \n\n");
	eslog::endln("BUILDER: MESH BUILT");
}

void buildOrderedFVM(InputMesh<OrderedUniqueNodes, OrderedUniqueFaces, OrderedRegions> &input, Mesh &mesh)
{
	eslog::startln("BUILDER: PROCESS FACED MESH", "BUILDER");
	mesh.dimension = 3;
	if (info::mpi::size == 1) {
		TemporalSequentialMesh<MergedNodes, OrderedFacesBalanced> grouped;
		TemporalSequentialMesh<MergedNodes, MergedElements> prepared;

		initializeSequentialFVM(input, grouped);
		eslog::checkpointln("BUILDER: DATA INITIALIZED");

		buildElementsFromFaces(grouped, prepared);
		eslog::checkpointln("BUILDER: ELEMENTS CREATED");

		fillSequentialMesh(prepared, *input.regions, mesh);
	} else {
		TemporalMesh<OrderedNodesBalanced, OrderedFacesBalanced> grouped;
		TemporalMesh<OrderedNodesBalanced, OrderedElementsBalanced> ordered;
		TemporalMesh<MergedNodes, MergedElements> clustered;
		TemporalMesh<LinkedNodes, MergedElements> linked;

		balanceFVM(input, grouped);
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
