
#include "builder.h"
#include "builder.utils.h"

#include "basis/containers/serializededata.h"
#include "basis/sfc/hilbertcurve.h"
#include "basis/utilities/packing.h"
#include "basis/utilities/parser.h"
#include "basis/utilities/utils.h"
#include "esinfo/eslog.hpp"
#include "esinfo/envinfo.h"
#include "mesh/store/elementsregionstore.h"
#include "mesh/store/boundaryregionstore.h"
#include "mesh/store/elementstore.h"
#include "mesh/store/nodestore.h"
#include "wrappers/mpi/communication.h"

#include <numeric>
#include <algorithm>

namespace espreso {

struct ChunkedVariables {
	ivector<esint> nperm, eperm;
	std::vector<esint> ndistribution, edistribution;

	esint datasize, nvsize, evsize;
	std::vector<esint> variableRequests, variableRequestsMap;
};

struct VariableLoader {
	std::vector<ElementData*> edata;
	std::vector<NodeData*> ndata;

	ChunkedVariables chunked;
};

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

void buildDecomposedFEM(NodesDomain &nodes, Elements &elements, OrderedRegions &regions, Mesh &mesh)
{
	eslog::startln("BUILDER: PROCESS DECOMPOSED MESH", "BUILDER");
	eslog::param("Nodes", nodes.coordinates.size());
	eslog::param("Elements", elements.etype.size());
	eslog::param("Total[B]", size(nodes) + size(elements));
	eslog::ln();

	eslog::info(" ================================= DECOMPOSED FEM BUILDER ====================== %12.3f s\n", eslog::duration());
	eslog::info(" ============================================================================================= \n");

	LinkedNodes linked;
	MergedElements merged;

	// 1. -> 8.
	// remove dangling
	// sort nodes?
	trivialUpdate(nodes, linked);
	trivialUpdate(elements, merged);

	// 9.
	fillNodes(linked, regions, mesh);
	fillElements(merged, regions, mesh);

	regionInfo(mesh);
	eslog::info(" ============================================================================================= \n\n");
	eslog::endln("BUILDER: MESH BUILT");
}

void buildDecomposedFVM(NodesDomain &nodes, Faces &faces, OrderedRegions &regions, Mesh &mesh)
{
	eslog::startln("BUILDER: PROCESS DECOMPOSED MESH", "BUILDER");
	eslog::param("Nodes", nodes.coordinates.size());
	eslog::param("Faces", faces.ftype.size());
	eslog::param("Total[B]", size(nodes) + size(faces));
	eslog::ln();

	eslog::info(" ================================== DECOMPOSED FVM BUILDER ===================== %12.3f s\n", eslog::duration());
	eslog::info(" ============================================================================================= \n");

	// 0. -> 1.
	Elements elements;
	buildElementsFromFaces(faces, elements);
	eslog::checkpointln("BUILDER: ELEMENTS BUILT");

	LinkedNodes linked;
	MergedElements merged;

	// 1. -> 8.
	// remove dangling
	// sort nodes?
	trivialUpdate(nodes, linked);
	trivialUpdate(elements, merged);

	rotateNormalsOut(linked, merged);
	eslog::checkpointln("BUILDER: ELEMENTS NODES REORDERED");

	// 9.
	fillNodes(linked, regions, mesh);
	fillElements(merged, regions, mesh);

	regionInfo(mesh);
	eslog::info(" ============================================================================================= \n\n");
	eslog::endln("BUILDER: MESH BUILT");
}

void buildOrderedFVM(NodesBlocks &nodes, FacesBlocks &faces, OrderedRegions &regions, Mesh &mesh)
{
	eslog::startln("BUILDER: PROCESS FACED MESH", "BUILDER");
	eslog::info(" ==================================== ORDERED FVM BUILDER ===================== %12.3f s\n", eslog::duration());
	eslog::info(" ============================================================================================= \n");
	mesh.dimension = 3;
	eslog::info(" == MESH DIMENSION %72d == \n", mesh.dimension);

	NodesHolder nh(std::move(nodes));
	FaceHolder fh(std::move(faces));
	ElementsHolder eh;

	if (info::mpi::size == 1) {
		// 0. -> 1.
		trivialUpdate(fh.blocks, fh.chunked);
		dynamic_cast<OrderedDistribution&>(eh.chunked) = fh.chunked.edist;
		buildElementsFromFaces(fh.chunked, eh.chunked);
		eslog::checkpointln("BUILDER: ELEMENTS BUILT");

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

		// 7. -> 8.
		globalToLocal(eh.clustered, eh.merged, nh.linked);
		eslog::checkpointln("BUILDER: ELEMENTS NODES REINDEXED");

		rotateNormalsOut(nh.linked, eh.merged);
		eslog::checkpointln("BUILDER: ELEMENTS NODES REORDERED");
	} else {
		ivector<esint> nbuckets, ebuckets, splitters;

		// 0. -> 1.
		balanceFVM(nh.blocks, fh.blocks, nh.balanced, fh.balanced);
		eslog::checkpointln("BUILDER: DATA BALANCED");

		dynamic_cast<BalancedDistribution&>(eh.balanced) = fh.balanced.edist;
		buildElementsFromFaces(fh.balanced, eh.balanced);
		eslog::param("OrderedNodesBalanced", nh.balanced.coordinates.size());
		eslog::param("OrderedElementsBalanced", eh.balanced.etype.size());
		eslog::param("OrderedTotalChunked[B]", size(nh.balanced) + size(nh.balanced));
		eslog::ln();
		eslog::checkpointln("BUILDER: ELEMENTS BUILT");

		// 3. synchronize mesh dimension and compute SFC
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
		globalToLocal(eh.clustered, eh.merged, nh.linked);
		eslog::checkpointln("BUILDER: ELEMENTS NODES REINDEXED");
		eslog::param("MergedElements", eh.merged.etype.size());
		eslog::ln();

		rotateNormalsOut(nh.linked, eh.merged);
		eslog::checkpointln("BUILDER: ELEMENTS NODES REORDERED");
	}

	// 9.
	fillNodes(nh.linked, regions, mesh);
	fillElements(eh.merged, regions, mesh);

	regionInfo(mesh);
	eslog::info(" ============================================================================================= \n\n");
	eslog::endln("BUILDER: MESH BUILT");
}

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
	buildElementsFromFaces(fh.chunked, eh.chunked);
	dynamic_cast<OrderedDistribution&>(eh.chunked) = fh.chunked.edist;
	eslog::checkpointln("BUILDER: ELEMENTS BUILT");

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

		// 7. -> 8.
		globalToLocal(eh.clustered, eh.merged, nh.linked);
		eslog::checkpointln("BUILDER: ELEMENTS NODES REINDEXED");

		rotateNormalsOut(nh.linked, eh.merged);
		eslog::checkpointln("BUILDER: ELEMENTS NODES REORDERED");
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
		eslog::checkpointln("BUILDER: ELEMENTS NODES REINDEXED");
		eslog::param("MergedElements", eh.merged.etype.size());
		eslog::ln();

		rotateNormalsOut(nh.linked, eh.merged);
		eslog::checkpointln("BUILDER: ELEMENTS NODES REORDERED");
	}

	// 9.
	fillNodes(nh.linked, regions, mesh);
	fillElements(eh.merged, regions, mesh);

	regionInfo(mesh);
	eslog::info(" ============================================================================================= \n\n");
	eslog::endln("BUILDER: MESH BUILT");
}

void orderedValuesInit(VariablesBlocks &variables, Mesh &mesh)
{
	variables.loader = new VariableLoader();
	ChunkedVariables &chunked = variables.loader->chunked;

	for (size_t v = 0; v < variables.nodes.size(); ++v) {
		if (StringCompare::caseInsensitiveEq("coordinates", variables.nodes[v].name)) {
			variables.loader->ndata.push_back(mesh.nodes->appendData(3, NamedData::DataType::VECTOR, "TRANSLATION")); // temporary store for coordinates
		} else {
			switch (variables.nodes[v].dimension) {
			case 1: variables.loader->ndata.push_back(mesh.nodes->appendData(1, NamedData::DataType::SCALAR, variables.nodes[v].name)); break;
			case 3: variables.loader->ndata.push_back(mesh.nodes->appendData(3, NamedData::DataType::VECTOR, variables.nodes[v].name)); break;
			default: break;
			}
		}
	}
	for (size_t v = 0; v < variables.elements.size(); ++v) {
		switch (variables.elements[v].dimension) {
		case 1: variables.loader->edata.push_back(mesh.elements->appendData(1, NamedData::DataType::SCALAR, variables.elements[v].name)); break;
		case 3: variables.loader->edata.push_back(mesh.elements->appendData(3, NamedData::DataType::VECTOR, variables.elements[v].name)); break;
		default: break;
		}
	}

//	chunked.nperm.resize(mesh.nodes->size);
	chunked.eperm.resize(mesh.elements->distribution.process.size);
//	std::iota(chunked.nperm.begin(), chunked.nperm.end(), 0);
//	std::sort(chunked.nperm.begin(), chunked.nperm.end(), [&] (esint i, esint j) { return (mesh.nodes->inputOffset->begin() + i)->front() < (mesh.nodes->inputOffset->begin() + j)->front(); });
	std::iota(chunked.eperm.begin(), chunked.eperm.end(), 0);
	std::sort(chunked.eperm.begin(), chunked.eperm.end(), [&] (esint i, esint j) { return (mesh.elements->inputOffset->begin() + i)->front() < (mesh.elements->inputOffset->begin() + j)->front(); });

	// there is only on block
//	chunked.ndistribution = Communication::getDistribution(variables.ndist.blocks.back().size);
	chunked.edistribution = Communication::getDistribution<esint>(variables.elements.front().data.size());

	std::vector<esint> sBuffer;
	sBuffer.reserve(info::mpi::size * 5 + chunked.nperm.size() + chunked.eperm.size());
	auto nit = chunked.nperm.begin(), eit = chunked.eperm.begin();
	for (int t = 0; t < info::mpi::size; ++t) {
		auto nbegin = nit, ebegin = eit;
		while (nit != chunked.nperm.end() && (mesh.nodes->inputOffset->begin() + *nit)->front() < chunked.ndistribution[t + 1]) { ++nit; }
		while (eit != chunked.eperm.end() && (mesh.elements->inputOffset->begin() + *eit)->front() < chunked.edistribution[t + 1]) { ++eit; }

		sBuffer.push_back(5 + (nit - nbegin) + (eit - ebegin)); // size
		sBuffer.push_back(t); // target
		sBuffer.push_back(info::mpi::rank); // source
		sBuffer.push_back(nit - nbegin); // nodes
		sBuffer.push_back(eit - ebegin); // elements
		for (auto it = nbegin; it != nit; ++it) {
			sBuffer.push_back((mesh.nodes->inputOffset->begin() + *it)->front());
		}
		for (auto it = ebegin; it != eit; ++it) {
			sBuffer.push_back((mesh.elements->inputOffset->begin() + *it)->front());
		}
	}

	if (!Communication::allToAllWithDataSizeAndTarget(sBuffer, chunked.variableRequests)) {
		eslog::error("cannot exchange output offsets.\n");
	}

	chunked.variableRequestsMap.resize(info::mpi::size);
	for (size_t i = 0; i < chunked.variableRequests.size(); i += chunked.variableRequests[i]) {
		chunked.variableRequestsMap[chunked.variableRequests[i + 2]] = i;
	}

	chunked.datasize = chunked.nvsize = chunked.evsize = 0;
	for (size_t v = 0; v < variables.nodes.size(); ++v) {
		chunked.nvsize += variables.nodes[v].dimension;
	}
	for (size_t v = 0; v < variables.elements.size(); ++v) {
		chunked.evsize += variables.elements[v].dimension;
	}
	for (size_t i = 0; i < chunked.variableRequests.size(); i += chunked.variableRequests[i]) {
		chunked.datasize += utils::reinterpret_size<esint, esfloat>(chunked.nvsize * chunked.variableRequests[i + 3]);
		chunked.datasize += utils::reinterpret_size<esint, esfloat>(chunked.evsize * chunked.variableRequests[i + 4]);
	}

	eslog::checkpointln("VARIABLES LOADER: PERMUTATION EXCHANGED");
}

void orderedValuesNext(VariablesBlocks &variables, Mesh &mesh)
{
	ChunkedVariables &chunked = variables.loader->chunked;

	std::vector<esint> sBuffer, rBuffer;
	sBuffer.reserve(info::mpi::size * 5 + chunked.datasize);
	for (int r = 0; r < info::mpi::size; ++r) {
		esint i = chunked.variableRequestsMap[r];
		esint nsize = chunked.variableRequests[i + 3];
		esint esize = chunked.variableRequests[i + 4];
		sBuffer.push_back(5 + utils::reinterpret_size<esint, esfloat>(chunked.nvsize * nsize + chunked.evsize * esize));
		sBuffer.push_back(r);
		sBuffer.push_back(info::mpi::rank);
		sBuffer.push_back(nsize);
		sBuffer.push_back(esize);
		size_t size = sBuffer.size(); // it is useless if we called correctly sBuffer.reserve(...)
		sBuffer.resize(sBuffer.size() + utils::reinterpret_size<esint, esfloat>(chunked.nvsize * nsize + chunked.evsize * esize));
		esfloat *c = reinterpret_cast<esfloat*>(sBuffer.data() + size);
		for (size_t v = 0; v < variables.nodes.size(); ++v) {
			if (variables.nodes[v].data.size()) {
				for (esint n = 0; n < nsize; ++n) {
					for (int d = 0; d < variables.nodes[v].dimension; ++d) {
						*c++ = variables.nodes[v].data[variables.nodes[v].dimension * (chunked.variableRequests[i + 5 + n] - chunked.ndistribution[info::mpi::rank]) + d];
					}
				}
			} else {
				c += nsize * variables.nodes[v].dimension;
			}
		}
		for (size_t v = 0; v < variables.elements.size(); ++v) {
			if (variables.elements[v].data.size()) {
				for (esint e = 0; e < esize; ++e) {
					for (int d = 0; d < variables.elements[v].dimension; ++d) {
						*c++ = variables.elements[v].data[variables.elements[v].dimension * (chunked.variableRequests[i + 5 + e + nsize] - chunked.edistribution[info::mpi::rank]) + d];
					}
				}
			} else {
				c += esize * variables.elements[v].dimension;
			}
		}
	}

	if (!Communication::allToAllWithDataSizeAndTarget(sBuffer, rBuffer)) {
		eslog::error("cannot exchange output data.\n");
	}
	sBuffer.clear();
	eslog::checkpointln("VARIABLES LOADER: VARIABLES EXCHANGED");

	std::vector<esint> rmap(info::mpi::size);
	for (size_t i = 0; i < rBuffer.size(); i += rBuffer[i]) {
		rmap[rBuffer[i + 2]] = i;
	}

	auto nit = chunked.nperm.begin(), eit = chunked.eperm.begin();
	for (int r = 0; r < info::mpi::size; ++r) {
		esint i = rmap[r];
		esint nsize = rBuffer[i + 3];
		esint esize = rBuffer[i + 4];
		esfloat *c = reinterpret_cast<esfloat*>(rBuffer.data() + i + 5);
		for (size_t v = 0; v < variables.nodes.size(); ++v) {
			for (esint n = 0; n < nsize; ++n) {
				for (int d = 0; d < variables.nodes[v].dimension; ++d) {
					variables.loader->ndata[v]->data[*(nit + n) * variables.nodes[v].dimension + d] = *c++;
				}
			}
		}
		nit += nsize;
		for (size_t v = 0; v < variables.elements.size(); ++v) {
			for (esint e = 0; e < esize; ++e) {
				for (int d = 0; d < variables.elements[v].dimension; ++d) {
					variables.loader->edata[v]->data[*(eit + e) * variables.elements[v].dimension + d] = *c++;
				}
			}
		}
		eit += esize;
	}
	eslog::checkpointln("VARIABLES LOADER: VARIABLES ASSIGNED");

	for (size_t v = 0; v < variables.nodes.size(); ++v) {
		if (StringCompare::caseInsensitiveEq("coordinates", variables.nodes[v].name)) {
			if (variables.nodes[v].data.size()) {
				if (mesh.nodes->originCoordinates == nullptr) {
					mesh.nodes->originCoordinates = mesh.nodes->coordinates;
					mesh.nodes->coordinates = new serializededata<esint, Point>(*mesh.nodes->originCoordinates);
				}
				for (size_t n = 0; n < mesh.nodes->coordinates->datatarray().size(); ++n) {
					mesh.nodes->coordinates->datatarray()[n].x = variables.loader->ndata[v]->data[3 * n + 0];
					mesh.nodes->coordinates->datatarray()[n].y = variables.loader->ndata[v]->data[3 * n + 1];
					mesh.nodes->coordinates->datatarray()[n].z = variables.loader->ndata[v]->data[3 * n + 2];
					variables.loader->ndata[v]->data[3 * n + 0] = mesh.nodes->coordinates->datatarray()[n].x - mesh.nodes->originCoordinates->datatarray()[n].x;
					variables.loader->ndata[v]->data[3 * n + 1] = mesh.nodes->coordinates->datatarray()[n].y - mesh.nodes->originCoordinates->datatarray()[n].y;
					variables.loader->ndata[v]->data[3 * n + 2] = mesh.nodes->coordinates->datatarray()[n].z - mesh.nodes->originCoordinates->datatarray()[n].z;
				}
				eslog::checkpointln("VARIABLES LOADER: COORDINATES UPDATED");
				mesh.updateMesh();
			} else {
				std::fill(variables.loader->ndata[v]->data.begin(), variables.loader->ndata[v]->data.end(), 0);
			}
		}
	}
}

void orderedValuesFinish(VariablesBlocks &variables)
{
	delete variables.loader;
}

void chunkedValuesInit(VariablesBlocks &variables, Mesh &mesh)
{
	variables.loader = new VariableLoader();
	ChunkedVariables &chunked = variables.loader->chunked;

	for (size_t v = 0; v < variables.nodes.size(); ++v) {
		if (StringCompare::caseInsensitiveEq("coordinates", variables.nodes[v].name)) {
			variables.loader->ndata.push_back(mesh.nodes->appendData(3, NamedData::DataType::VECTOR, "TRANSLATION")); // temporary store for coordinates
		} else {
			switch (variables.nodes[v].dimension) {
			case 1: variables.loader->ndata.push_back(mesh.nodes->appendData(1, NamedData::DataType::SCALAR, variables.nodes[v].name)); break;
			case 3: variables.loader->ndata.push_back(mesh.nodes->appendData(3, NamedData::DataType::VECTOR, variables.nodes[v].name)); break;
			default: break;
			}
		}
		variables.loader->ndata.back()->filter = variables.nodes[v].filter;
	}
	for (size_t v = 0; v < variables.elements.size(); ++v) {
		switch (variables.elements[v].dimension) {
		case 1: variables.loader->edata.push_back(mesh.elements->appendData(1, NamedData::DataType::SCALAR, variables.elements[v].name)); break;
		case 3: variables.loader->edata.push_back(mesh.elements->appendData(3, NamedData::DataType::VECTOR, variables.elements[v].name)); break;
		default: break;
		}
		variables.loader->edata.back()->filter = variables.elements[v].filter;
	}

	chunked.nperm.resize(mesh.nodes->size);
	chunked.eperm.resize(mesh.elements->distribution.process.size);
	std::iota(chunked.nperm.begin(), chunked.nperm.end(), 0);
	std::sort(chunked.nperm.begin(), chunked.nperm.end(), [&] (esint i, esint j) { return (mesh.nodes->inputOffset->begin() + i)->front() < (mesh.nodes->inputOffset->begin() + j)->front(); });
	std::iota(chunked.eperm.begin(), chunked.eperm.end(), 0);
	std::sort(chunked.eperm.begin(), chunked.eperm.end(), [&] (esint i, esint j) { return (mesh.elements->inputOffset->begin() + i)->front() < (mesh.elements->inputOffset->begin() + j)->front(); });

	// there is only on block
	chunked.ndistribution = Communication::getDistribution(variables.ndist.blocks.back().size);
	chunked.edistribution = Communication::getDistribution(variables.edist.blocks.back().size);

	std::vector<esint> sBuffer;
	sBuffer.reserve(info::mpi::size * 5 + chunked.nperm.size() + chunked.eperm.size());
	auto nit = chunked.nperm.begin(), eit = chunked.eperm.begin();
	for (int t = 0; t < info::mpi::size; ++t) {
		auto nbegin = nit, ebegin = eit;
		while (nit != chunked.nperm.end() && (mesh.nodes->inputOffset->begin() + *nit)->front() < chunked.ndistribution[t + 1]) { ++nit; }
		while (eit != chunked.eperm.end() && (mesh.elements->inputOffset->begin() + *eit)->front() < chunked.edistribution[t + 1]) { ++eit; }

		sBuffer.push_back(5 + (nit - nbegin) + (eit - ebegin)); // size
		sBuffer.push_back(t); // target
		sBuffer.push_back(info::mpi::rank); // source
		sBuffer.push_back(nit - nbegin); // nodes
		sBuffer.push_back(eit - ebegin); // elements
		for (auto it = nbegin; it != nit; ++it) {
			sBuffer.push_back((mesh.nodes->inputOffset->begin() + *it)->front());
		}
		for (auto it = ebegin; it != eit; ++it) {
			sBuffer.push_back((mesh.elements->inputOffset->begin() + *it)->front());
		}
	}

	if (!Communication::allToAllWithDataSizeAndTarget(sBuffer, chunked.variableRequests)) {
		eslog::error("cannot exchange output offsets.\n");
	}

	chunked.variableRequestsMap.resize(info::mpi::size);
	for (size_t i = 0; i < chunked.variableRequests.size(); i += chunked.variableRequests[i]) {
		chunked.variableRequestsMap[chunked.variableRequests[i + 2]] = i;
	}

	chunked.datasize = chunked.nvsize = chunked.evsize = 0;
	for (size_t v = 0; v < variables.nodes.size(); ++v) {
		chunked.nvsize += variables.nodes[v].dimension;
	}
	for (size_t v = 0; v < variables.elements.size(); ++v) {
		chunked.evsize += variables.elements[v].dimension;
	}
	for (size_t i = 0; i < chunked.variableRequests.size(); i += chunked.variableRequests[i]) {
		chunked.datasize += utils::reinterpret_size<esint, esfloat>(chunked.nvsize * chunked.variableRequests[i + 3]);
		chunked.datasize += utils::reinterpret_size<esint, esfloat>(chunked.evsize * chunked.variableRequests[i + 4]);
	}

	eslog::checkpointln("VARIABLES LOADER: PERMUTATION EXCHANGED");
}

void chunkedValuesNext(VariablesBlocks &variables, Mesh &mesh)
{
	ChunkedVariables &chunked = variables.loader->chunked;

	std::vector<esint> sBuffer, rBuffer;
	sBuffer.reserve(info::mpi::size * 5 + chunked.datasize);
	for (int r = 0; r < info::mpi::size; ++r) {
		esint i = chunked.variableRequestsMap[r];
		esint nsize = chunked.variableRequests[i + 3];
		esint esize = chunked.variableRequests[i + 4];
		sBuffer.push_back(5 + utils::reinterpret_size<esint, esfloat>(chunked.nvsize * nsize + chunked.evsize * esize));
		sBuffer.push_back(r);
		sBuffer.push_back(info::mpi::rank);
		sBuffer.push_back(nsize);
		sBuffer.push_back(esize);
		size_t size = sBuffer.size(); // it is useless if we called correctly sBuffer.reserve(...)
		sBuffer.resize(sBuffer.size() + utils::reinterpret_size<esint, esfloat>(chunked.nvsize * nsize + chunked.evsize * esize));
		esfloat *c = reinterpret_cast<esfloat*>(sBuffer.data() + size);
		for (size_t v = 0; v < variables.nodes.size(); ++v) {
			if (variables.nodes[v].data.size()) {
				for (esint n = 0; n < nsize; ++n) {
					for (int d = 0; d < variables.nodes[v].dimension; ++d) {
						*c++ = variables.nodes[v].data[variables.nodes[v].dimension * (chunked.variableRequests[i + 5 + n] - chunked.ndistribution[info::mpi::rank]) + d];
					}
				}
			} else {
				c += nsize * variables.nodes[v].dimension;
			}
		}
		for (size_t v = 0; v < variables.elements.size(); ++v) {
			if (variables.elements[v].data.size()) {
				for (esint e = 0; e < esize; ++e) {
					for (int d = 0; d < variables.elements[v].dimension; ++d) {
						*c++ = variables.elements[v].data[variables.elements[v].dimension * (chunked.variableRequests[i + 5 + e + nsize] - chunked.edistribution[info::mpi::rank]) + d];
					}
				}
			} else {
				c += esize * variables.elements[v].dimension;
			}
		}
	}

	if (!Communication::allToAllWithDataSizeAndTarget(sBuffer, rBuffer)) {
		eslog::error("cannot exchange output data.\n");
	}
	sBuffer.clear();
	eslog::checkpointln("VARIABLES LOADER: VARIABLES EXCHANGED");

	std::vector<esint> rmap(info::mpi::size);
	for (size_t i = 0; i < rBuffer.size(); i += rBuffer[i]) {
		rmap[rBuffer[i + 2]] = i;
	}

	auto nit = chunked.nperm.begin(), eit = chunked.eperm.begin();
	for (int r = 0; r < info::mpi::size; ++r) {
		esint i = rmap[r];
		esint nsize = rBuffer[i + 3];
		esint esize = rBuffer[i + 4];
		esfloat *c = reinterpret_cast<esfloat*>(rBuffer.data() + i + 5);
		for (size_t v = 0; v < variables.nodes.size(); ++v) {
			for (esint n = 0; n < nsize; ++n) {
				for (int d = 0; d < variables.nodes[v].dimension; ++d) {
					variables.loader->ndata[v]->data[*(nit + n) * variables.nodes[v].dimension + d] = *c++;
				}
			}
		}
		nit += nsize;
		for (size_t v = 0; v < variables.elements.size(); ++v) {
			for (esint e = 0; e < esize; ++e) {
				for (int d = 0; d < variables.elements[v].dimension; ++d) {
					variables.loader->edata[v]->data[*(eit + e) * variables.elements[v].dimension + d] = *c++;
				}
			}
		}
		eit += esize;
	}
	eslog::checkpointln("VARIABLES LOADER: VARIABLES ASSIGNED");

	for (size_t v = 0; v < variables.nodes.size(); ++v) {
		if (StringCompare::caseInsensitiveEq("coordinates", variables.nodes[v].name)) {
			if (variables.nodes[v].data.size()) {
				if (mesh.nodes->originCoordinates == nullptr) {
					mesh.nodes->originCoordinates = mesh.nodes->coordinates;
					mesh.nodes->coordinates = new serializededata<esint, Point>(*mesh.nodes->originCoordinates);
				}
				for (size_t n = 0; n < mesh.nodes->coordinates->datatarray().size(); ++n) {
					mesh.nodes->coordinates->datatarray()[n].x = variables.loader->ndata[v]->data[3 * n + 0];
					mesh.nodes->coordinates->datatarray()[n].y = variables.loader->ndata[v]->data[3 * n + 1];
					mesh.nodes->coordinates->datatarray()[n].z = variables.loader->ndata[v]->data[3 * n + 2];
					variables.loader->ndata[v]->data[3 * n + 0] = mesh.nodes->coordinates->datatarray()[n].x - mesh.nodes->originCoordinates->datatarray()[n].x;
					variables.loader->ndata[v]->data[3 * n + 1] = mesh.nodes->coordinates->datatarray()[n].y - mesh.nodes->originCoordinates->datatarray()[n].y;
					variables.loader->ndata[v]->data[3 * n + 2] = mesh.nodes->coordinates->datatarray()[n].z - mesh.nodes->originCoordinates->datatarray()[n].z;
				}
				eslog::checkpointln("VARIABLES LOADER: COORDINATES UPDATED");
				mesh.updateMesh();
			} else {
				std::fill(variables.loader->ndata[v]->data.begin(), variables.loader->ndata[v]->data.end(), 0);
			}
		}
	}
}

void chunkedValuesFinish(VariablesBlocks &variables)
{
	delete variables.loader;
}

}
}
