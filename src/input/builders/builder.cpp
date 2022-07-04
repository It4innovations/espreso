
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

int getDimension(OrderedElementsBalanced &balanced)
{
	int dimension = 0;
	for (esint e = 0; e < balanced.size && dimension < 3; ++e) {
		switch (Mesh::element(balanced.etype[e]).type) {
		case Element::TYPE::POINT:  dimension = std::max(0, dimension); break;
		case Element::TYPE::LINE:   dimension = std::max(1, dimension); break;
		case Element::TYPE::PLANE:  dimension = std::max(2, dimension); break;
		case Element::TYPE::VOLUME: dimension = std::max(3, dimension); break;
		}
	}
	Communication::allReduce(&dimension, nullptr, 1, MPI_INT, MPI_MAX);
	eslog::info(" == MESH DIMENSION %72d == \n", dimension);
	return dimension;
}

static void regionInfo(Mesh &mesh)
{
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
}

void buildOrderedFEM(OrderedNodes &nodes, OrderedElements &elements, OrderedRegions &regions, Mesh &mesh)
{
	eslog::startln("BUILDER: PROCESS ORDERED MESH", "BUILDER");
	eslog::param("OrderedNodes", size(nodes));
	eslog::param("OrderedElements", size(elements));
	eslog::param("OrderedTotal", size(nodes) + size(elements));

	eslog::info(" ==================================== ORDERED FEM BUILDER ===================== %12.3f s\n", eslog::duration());
	eslog::info(" ============================================================================================= \n");

	NodesHolder nh(std::move(nodes));
	ElementsHolder eh(std::move(elements));

	/**
	 * Standard workflow for ordered database with regions.
	 *
	 *    +-- NODES --+- ELEMENS -+
	 *    |-----------+-----------| <-- parallel parser provides data
	 * 1. | ORDERED   | ORDERED   | <-- initial setting for the builder
	 * 2. | BALANCED  | BALANCED  |
	 * 3. |-----------------------| <-- compute SFC and assign buckets
	 * 4. | CLUSTERED | CLUSTERED |
	 * 5. |-----------------------| <-- approximate neighbors
	 * 6. | MERGED    | CLUSTERED |
	 * 7. | LINKED    | CLUSTERED |
	 * 8. | LINKED    | MERGED    |
	 * 9. +-----------------------+ <-- fill the mesh
	 */
	if (info::mpi::size == 1) {
		// 1. -> 2. trivial
		trivialUpdate(nh.ordered, nh.balanced);
		trivialUpdate(eh.ordered, eh.balanced);
		mesh.dimension = getDimension(eh.balanced);

		// 2. -> 4. trivial
		trivialUpdate(nh.balanced, nh.clustered);
		trivialUpdate(eh.balanced, eh.clustered);
		eslog::checkpointln("BUILDER: DATA INITIALIZED");

		eslog::info(" ==    -   -   -   -   -   -   -   -   -   -   -   -   -   -   -   -   -   -   -   -   -    == \n");
		// 4. skipped
		// 5. -> 6. just local search
		searchDuplicatedNodes(nh.clustered, nh.merged);
		eslog::checkpointln("BUILDER: DUPLICATED NODES FOUND");

		// 6. -> 7. trivial
		mergeDuplicatedNodes(nh.merged); // fix 'linkup' side effect
		trivialUpdate(nh.merged, nh.linked);

		// 7. -> 8. just local search
		mergeDuplicatedElements(eh.clustered, eh.merged, nh.linked, mesh.dimension);
		eslog::checkpointln("BUILDER: DUPLICATED ELEMENTS FOUND");
	} else {
		ivector<esint> nbuckets, ebuckets, splitters;

		// 1. -> 2.
		balanceFEM(nh.ordered, eh.ordered, nh.balanced, eh.balanced);
		eslog::checkpointln("BUILDER: DATA BALANCED");
		eslog::param("OrderedNodesBalanced", size(nh.balanced));
		eslog::param("OrderedElementsBalanced", size(eh.balanced));
		eslog::param("OrderedTotalBalanced", size(nh.balanced) + size(nh.balanced));

		// 3. synchronize mesh dimension and compute SFC
		mesh.dimension = getDimension(eh.balanced);
		HilbertCurve<esfloat> sfc(mesh.dimension, SFCDEPTH, nh.balanced.coordinates.size(), nh.balanced.coordinates.data());
		eslog::checkpointln("BUILDER: SFC BUILT");

		assignBuckets(nh.balanced, eh.balanced, sfc, nbuckets, ebuckets);
		eslog::checkpointln("BUILDER: SFC BUCKETS ASSIGNED");

		// 2. -> 4.
		clusterize(nh.balanced, eh.balanced, nbuckets, ebuckets, sfc.buckets(sfc.depth), nh.clustered, eh.clustered, splitters);
		utils::clearVectors(nbuckets, ebuckets);
		eslog::checkpointln("BUILDER: MESH CLUSTERIZED");
		eslog::param("ClusteredNodes", size(nh.clustered));
		eslog::param("ClusteredElements", size(eh.clustered));

		// 5.
		computeSFCNeighbors(sfc, splitters, nh.linked.neighbors); // neighbors are approximated here
		eslog::checkpointln("BUILDER: NEIGHBORS APPROXIMATED");
		eslog::info(" ==    -   -   -   -   -   -   -   -   -   -   -   -   -   -   -   -   -   -   -   -   -    == \n");

		// 4. -> 6.
		exchangeSFCBoundaryNodes(sfc, splitters, nh.linked.neighbors, nh.clustered);
		searchDuplicatedNodes(nh.clustered, nh.merged);
		utils::clearVector(splitters);
		eslog::checkpointln("BUILDER: DUPLICATED NODES FOUND");

		// 6. -> 7.
		linkup(nh.merged, nh.linked, eh.clustered);
		eslog::checkpointln("BUILDER: LINKED UP");
		eslog::param("LinkedNodes", size(nh.linked));

		// 7. -> 8.
		mergeDuplicatedElements(eh.clustered, eh.merged, nh.linked, mesh.dimension);
		eslog::checkpointln("BUILDER: DUPLICATED ELEMENTS FOUND");
		eslog::param("MergedElements", size(eh.merged));
	}
	// 9.
	fillNodes(nh.linked, regions, mesh);
	fillElements(eh.merged, regions, mesh);

	regionInfo(mesh);
	eslog::info(" ============================================================================================= \n\n");
	eslog::endln("BUILDER: MESH BUILT");
}

//void buildOrderedFVM(InputMesh<OrderedNodes, OrderedFaces, OrderedRegions> &input, Mesh &mesh)
//{
//	eslog::startln("BUILDER: PROCESS FACED MESH", "BUILDER");
//	mesh.dimension = 3;
//	if (info::mpi::size == 1) {
//		TemporalSequentialMesh<MergedNodes, OrderedFacesBalanced> grouped;
//		TemporalSequentialMesh<MergedNodes, MergedElements> prepared;
//
//		initializeSequentialFVM(input, grouped);
//		eslog::checkpointln("BUILDER: DATA INITIALIZED");
//
//		buildElementsFromFaces(grouped, prepared);
//		eslog::checkpointln("BUILDER: ELEMENTS CREATED");
//
//		fillSequentialMesh(prepared, *input.regions, mesh);
//	} else {
//		TemporalMesh<OrderedNodesBalanced, OrderedFacesBalanced> grouped;
//		TemporalMesh<OrderedNodesBalanced, OrderedElementsBalanced> ordered;
//		TemporalMesh<MergedNodes, MergedElements> clustered;
//		TemporalMesh<LinkedNodes, MergedElements> linked;
//
//		balanceFVM(input, grouped);
//		eslog::checkpointln("BUILDER: DATA BALANCED");
//
//		buildElementsFromFaces(grouped, ordered);
//		eslog::checkpointln("BUILDER: ELEMENTS CREATED");
//
//		ivector<esint> nbuckets, ebuckets, splitters;
//		HilbertCurve<esfloat> sfc(mesh.dimension, SFCDEPTH, ordered.nodes->coordinates.size(), ordered.nodes->coordinates.data());
//		eslog::checkpointln("BUILDER: SFC BUILT");
//
//		assignBuckets(ordered, sfc, nbuckets, ebuckets);
//		eslog::checkpointln("BUILDER: SFC BUCKETS ASSIGNED");
//
//		clusterize(ordered, nbuckets, ebuckets, sfc.buckets(sfc.depth), clustered, splitters);
//		utils::clearVectors(nbuckets, ebuckets);
//		eslog::checkpointln("BUILDER: MESH CLUSTERIZED");
//
//		computeSFCNeighbors(sfc, splitters, linked.nodes->neighbors); // neighbors are approximated here
//		eslog::checkpointln("BUILDER: NEIGHBORS APPROXIMATED");
//
//		linkup(clustered, linked);
//		reindexToLocal(linked);
//		eslog::checkpointln("BUILDER: LINKED UP");
//
//		fillMesh(linked, *input.regions, mesh);
//	}
//	eslog::endln("BUILDER: MESH BUILT");
//}

void buildDecomposedFVM(OrderedNodes &nodes, OrderedFaces &elements, OrderedRegions &regions, Mesh &mesh)
{
	eslog::startln("BUILDER: PROCESS FACED MESH", "BUILDER");
	eslog::info(" ================================== DECOMPOSED FVM BUILDER ==================== %12.3f s\n", eslog::duration());
	eslog::info(" ============================================================================================= \n");
	mesh.dimension = 3;
	eslog::info(" == MESH DIMENSION %72d == \n", mesh.dimension);



	if (info::mpi::size == 1) {

	}

	regionInfo(mesh);
	eslog::info(" ============================================================================================= \n\n");
	eslog::endln("BUILDER: MESH BUILT");
}
//
//void buildDecomposedFVM(InputMesh<OrderedNodes, OrderedFaces, OrderedRegions> &input, Mesh &mesh)
//{
//	eslog::startln("BUILDER: PROCESS FACED MESH", "BUILDER");
//	eslog::info(" ================================== DECOMPOSED FVM BUILDER ==================== %12.3f s\n", eslog::duration());
//	eslog::info(" ============================================================================================= \n");
//	mesh.dimension = 3;
//	eslog::info(" == MESH DIMENSION %72d == \n", mesh.dimension);
//
//	if (info::mpi::size == 1) {
//		TemporalSequentialMesh<ClusteredNodes, OrderedFacesBalanced> grouped;
//		TemporalSequentialMesh<MergedNodes, OrderedFacesBalanced> merged;
//		TemporalSequentialMesh<MergedNodes, ClusteredElements> fem;
//		TemporalSequentialMesh<MergedNodes, MergedElements> prepared;
//
//		initializeSequentialFVM(input, grouped);
//		eslog::checkpointln("BUILDER: DATA INITIALIZED");
//
//		searchDuplicatedNodes(grouped, merged);
//		eslog::info(" == DUPLICATED NODES %70d == \n", merged.nodes->duplication.size());
//		eslog::checkpointln("BUILDER: DUPLICATED NODES FOUND");
//
//		buildElementsFromFaces(merged, fem);
//		eslog::checkpointln("BUILDER: ELEMENTS CREATED");
//
//		searchDuplicatedElements(fem, prepared, mesh.dimension);
//		eslog::info(" == DUPLICATED ELEMENTS %67d == \n", prepared.elements->duplication.size());
//		eslog::checkpointln("BUILDER: DUPLICATED ELEMENTS FOUND");
//
//		fillSequentialMesh(prepared, *input.regions, mesh);
//	} else {
//		TemporalMesh<OrderedNodesBalanced, OrderedFacesBalanced> grouped;
//		TemporalMesh<OrderedNodesBalanced, OrderedElementsBalanced> ordered;
//		TemporalMesh<MergedNodes, MergedElements> clustered;
//		TemporalMesh<LinkedNodes, MergedElements> linked;

//		balanceFVM(input, grouped);
//		eslog::checkpointln("BUILDER: DATA BALANCED");
//
//		buildElementsFromFaces(grouped, ordered);
//		eslog::checkpointln("BUILDER: ELEMENTS CREATED");
//
//		ivector<esint> nbuckets, ebuckets, splitters;
//		HilbertCurve<esfloat> sfc(mesh.dimension, SFCDEPTH, ordered.nodes->coordinates.size(), ordered.nodes->coordinates.data());
//		eslog::checkpointln("BUILDER: SFC BUILT");
//
//		assignBuckets(ordered, sfc, nbuckets, ebuckets);
//		eslog::checkpointln("BUILDER: SFC BUCKETS ASSIGNED");
//
//		clusterize(ordered, nbuckets, ebuckets, sfc.buckets(sfc.depth), clustered, splitters);
//		utils::clearVectors(nbuckets, ebuckets);
//		eslog::checkpointln("BUILDER: MESH CLUSTERIZED");
//
//		computeSFCNeighbors(sfc, splitters, linked.nodes->neighbors); // neighbors are approximated here
//		eslog::checkpointln("BUILDER: NEIGHBORS APPROXIMATED");
//
//		linkup(clustered, linked);
//		reindexToLocal(linked);
//		eslog::checkpointln("BUILDER: LINKED UP");
//
//		fillMesh(linked, *input.regions, mesh);
//	}
//
//	regionInfo(mesh);
//	eslog::info(" ============================================================================================= \n\n");
//	eslog::endln("BUILDER: MESH BUILT");
//}

}
}
