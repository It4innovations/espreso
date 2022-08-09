
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

int getDimension(Elements &elements)
{
	int dimension = 0;
	for (size_t e = 0; e < elements.etype.size() && dimension < 3; ++e) {
		switch (Mesh::element(elements.etype[e]).type) {
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

void buildOrderedFEM(NodesBlocks &nodes, ElementsBlocks &elements, OrderedRegions &regions, Mesh &mesh)
{
	eslog::startln("BUILDER: PROCESS ORDERED MESH", "BUILDER");
	eslog::param("OrderedNodes", nodes.coordinates.size());
	eslog::param("OrderedElements", elements.etype.size());
	eslog::param("OrderedTotal[B]", size(nodes) + size(elements));
	eslog::ln();

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
		trivialUpdate(nh.blocks, nh.chunked);
		trivialUpdate(eh.blocks, eh.chunked);
		mesh.dimension = getDimension(eh.chunked);

		// 2. -> 4. trivial
		trivialUpdate(nh.chunked, nh.clustered);
		trivialUpdate(eh.chunked, eh.clustered);
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
		balanceFEM(nh.blocks, eh.blocks, nh.balanced, eh.balanced);
		eslog::checkpointln("BUILDER: DATA BALANCED");
		eslog::param("OrderedNodesBalanced", nh.balanced.coordinates.size());
		eslog::param("OrderedElementsBalanced", eh.balanced.etype.size());
		eslog::param("OrderedTotalBalanced[B]", size(nh.balanced) + size(nh.balanced));
		eslog::ln();

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
		eslog::param("ClusteredNodes", nh.clustered.coordinates.size());
		eslog::param("ClusteredElements", eh.clustered.etype.size());
		eslog::ln();

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
		eslog::param("LinkedNodes", nh.linked.coordinates.size());
		eslog::ln();

		// 7. -> 8.
		mergeDuplicatedElements(eh.clustered, eh.merged, nh.linked, mesh.dimension);
		eslog::checkpointln("BUILDER: DUPLICATED ELEMENTS FOUND");
		eslog::param("MergedElements", eh.merged.etype.size());
		eslog::ln();
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

void buildChunkedFVM(NodesBlocks &nodes, FacesBlocks &faces, OrderedRegions &regions, Mesh &mesh)
{
	eslog::startln("BUILDER: PROCESS FACED MESH", "BUILDER");
	eslog::info(" ================================== DECOMPOSED FVM BUILDER ==================== %12.3f s\n", eslog::duration());
	eslog::info(" ============================================================================================= \n");
	mesh.dimension = 3;
	eslog::info(" == MESH DIMENSION %72d == \n", mesh.dimension);

	NodesHolder nh(std::move(nodes));
	FaceHolder fh(std::move(faces));
	ElementsHolder eh;

	// 0. -> 1.
	trivialUpdate(fh.blocks, fh.chunked);
	buildElementsFromFaces(fh.chunked, eh.chunked, nh.blocks);

	if (info::mpi::size == 1) {
		// 1. -> 2. trivial
		trivialUpdate(nh.blocks, nh.chunked);

		// 2. -> 4. trivial
		trivialUpdate(nh.chunked, nh.clustered);
		trivialUpdate(eh.chunked, eh.clustered);
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
		globalToLocal(eh.clustered, eh.merged, nh.linked);
		eslog::checkpointln("BUILDER: DUPLICATED ELEMENTS FOUND");
	} else {
		ivector<esint> nbuckets, ebuckets, splitters;

		// 1. -> 2. trivial
		trivialUpdate(nh.blocks, nh.chunked);
		eslog::param("OrderedNodesChunked", nh.chunked.coordinates.size());
		eslog::param("OrderedElementsChunked", eh.chunked.etype.size());
		eslog::param("OrderedTotalChunked[B]", size(nh.chunked) + size(nh.chunked));
		eslog::ln();

		// 3. synchronize mesh dimension and compute SFC
		HilbertCurve<esfloat> sfc(mesh.dimension, SFCDEPTH, nh.chunked.coordinates.size(), nh.chunked.coordinates.data());
		eslog::checkpointln("BUILDER: SFC BUILT");

		assignBuckets(nh.chunked, eh.chunked, sfc, nbuckets, ebuckets);
		eslog::checkpointln("BUILDER: SFC BUCKETS ASSIGNED");

		// 2. -> 4.
		clusterize(nh.chunked, eh.chunked, nbuckets, ebuckets, sfc.buckets(sfc.depth), nh.clustered, eh.clustered, splitters);
		utils::clearVectors(nbuckets, ebuckets);
		eslog::checkpointln("BUILDER: MESH CLUSTERIZED");
		eslog::param("ClusteredNodes", nh.clustered.coordinates.size());
		eslog::param("ClusteredElements", eh.clustered.etype.size());
		eslog::ln();

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
		eslog::param("LinkedNodes", nh.linked.coordinates.size());
		eslog::ln();

		// 7. -> 8.
		globalToLocal(eh.clustered, eh.merged, nh.linked);
		eslog::checkpointln("BUILDER: DUPLICATED ELEMENTS FOUND");
		eslog::param("MergedElements", eh.merged.etype.size());
		eslog::ln();
	}

	// 9.
	fillNodes(nh.linked, regions, mesh);
	fillElements(eh.merged, regions, mesh);

	regionInfo(mesh);
	eslog::info(" ============================================================================================= \n\n");
	eslog::endln("BUILDER: MESH BUILT");
}

}
}
