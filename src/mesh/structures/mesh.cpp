
#include "mesh.h"
#include "mkl.h"
#include "cilk/cilk.h"

namespace espreso {

Mesh::Mesh():_elements(0)
{
	_partPtrs.resize(2);
	_partPtrs[0] = 0;
	_partPtrs[1] = 0;
}

void Mesh::partitiate(size_t parts)
{
	if (parts == 1 && this->parts() == 1) {
		_partPtrs = { 0, (eslocal)_elements.size() };
		mapElementsToDomains();
	} else {
		std::vector<eslocal> ePartition = getPartition(0, _elements.size(), parts);

		_partPtrs = std::vector<eslocal>(parts + 1, 0);
		for (size_t i = 0; i < _elements.size(); i++) {
			_elements[i]->domains().clear();
			_elements[i]->domains().push_back(ePartition[i]);
			_partPtrs[ePartition[i]]++;
		}

		std::sort(_elements.begin(), _elements.end(), [] (const Element* e1, const Element* e2) { return e1->domains()[0] < e2->domains()[0]; });
		ESTEST(MANDATORY) << "subdomain without element" << (std::any_of(_partPtrs.begin(), _partPtrs.end() - 1, [] (eslocal size) { return size == 0; }) ? TEST_FAILED : TEST_PASSED);
		Esutils::sizesToOffsets(_partPtrs);
	}

	_fixPoints.clear();
	_corners.clear();
	mapFacesToDomains();
	mapEdgesToDomains();
	mapNodesToDomains();
	mapCoordinatesToDomains();
}

void APIMesh::partitiate(size_t parts)
{
	if (parts == 1 && this->parts() == 1) {
		_partPtrs = { 0, (eslocal)_elements.size() };
		mapElementsToDomains();
	} else {
		std::vector<eslocal> ePartition = getPartition(0, _elements.size(), parts);

		_partPtrs = std::vector<eslocal>(parts + 1, 0);
		for (size_t i = 0; i < _elements.size(); i++) {
			_elements[i]->domains().clear();
			_elements[i]->domains().push_back(ePartition[i]);
			_partPtrs[ePartition[i]]++;
		}

		std::iota(_permutation.begin(), _permutation.end(), 0);
		std::sort(_permutation.begin(), _permutation.end(), [&] (eslocal i, eslocal j) { return _elements[i]->domains()[0] < _elements[j]->domains()[0]; });

		std::sort(_elements.begin(), _elements.end(), [] (const Element* e1, const Element* e2) { return e1->domains()[0] < e2->domains()[0]; });
		ESTEST(MANDATORY) << "subdomain without element" << (std::any_of(_partPtrs.begin(), _partPtrs.end() - 1, [] (eslocal size) { return size == 0; }) ? TEST_FAILED : TEST_PASSED);
		Esutils::sizesToOffsets(_partPtrs);
	}

	mapNodesToDomains();
	mapDOFsToDomains();
	mapCoordinatesToDomains();
}

void Mesh::computeFixPoints(size_t number)
{
	if (_fixPoints.size() && _fixPoints[0].size() == number) {
		return;
	}

	_fixPoints.resize(parts());

	cilk_for (size_t part = 0; part < parts(); part++) {
		size_t max = _partPtrs[part + 1] - _partPtrs[part];
		size_t points = std::min(number, max);
		std::vector<eslocal> fixPoints(points);
		std::vector<eslocal> eSubPartition = getPartition(_partPtrs[part], _partPtrs[part + 1], points);

		for (eslocal j = 0; j < points; j++) {
			fixPoints[j] = getCentralNode(_partPtrs[part], _partPtrs[part + 1], eSubPartition, part, j);
		}
		std::sort(fixPoints.begin(), fixPoints.end());

		// Remove the same points
		Esutils::removeDuplicity(fixPoints);

		_fixPoints[part].reserve(fixPoints.size());
		for (size_t i = 0; i < fixPoints.size(); i++) {
			_fixPoints[part].push_back(_nodes[fixPoints[i]]);
		}
	}
}

static void checkMETISResult(eslocal result)
{
	switch (result) {
	case METIS_ERROR_INPUT:
		ESINFO(ERROR) << "An input for METIS procedure is incorrect.\n";
	case METIS_ERROR_MEMORY:
		ESINFO(ERROR) << "There is not enough memory for compute a partition.\n";
	case METIS_ERROR:
		ESINFO(ERROR) << "METIS fail computation.\n";
	}
}

static std::vector<eslocal> computePartition(std::vector<eslocal> &elements, std::vector<eslocal> &nodes, eslocal eSize, eslocal nSize, eslocal nCommon, eslocal parts)
{
	if (parts == 1) {
		return std::vector<eslocal> (elements.size(), 0);
	}
	// INPUTS
	eslocal options[METIS_NOPTIONS];
	METIS_SetDefaultOptions(options);
	options[METIS_OPTION_CONTIG] = 1;

	// OUTPUTS
	eslocal objval;
	std::vector<eslocal> ePartition(eSize), nPartition(nSize);

	checkMETISResult(METIS_PartMeshDual(
			&eSize,
			&nSize,
			elements.data(),
			nodes.data(),
			NULL,		// weights of nodes
			NULL,		// size of nodes
			&nCommon,
			&parts,
			NULL,		// weights of parts
			options,
			&objval,
			ePartition.data(),
			nPartition.data()));

	return ePartition;
}

static eslocal localIndex(const std::vector<eslocal> &nodeProjection, eslocal index)
{
	return std::lower_bound(nodeProjection.begin(), nodeProjection.end(), index) - nodeProjection.begin();
}

static void METISMeshRepresentation(
		const std::vector<Element*> &elements, size_t begin, size_t end,
		std::vector<eslocal> &nodeProjection, std::vector<eslocal> &metisElements, std::vector<eslocal> &metisNodes, eslocal &nCommon)
{
	for (size_t e = begin; e < end; e++) {
		for (size_t n = 0; n < elements[e]->coarseNodes(); n++) {
			nodeProjection.push_back(elements[e]->node(n));
		}
	}

	std::sort(nodeProjection.begin(), nodeProjection.end());
	Esutils::removeDuplicity(nodeProjection);

	// create array storing pointers to elements' nodes
	metisElements.reserve(end - begin + 1);
	metisElements.push_back(0);
	for (eslocal i = begin; i < end; i++) {
		metisElements.push_back(metisElements.back() + elements[i]->coarseNodes());
	}

	// create array of nodes
	metisNodes.reserve(metisElements.back());
	// number of common nodes to be neighbor
	nCommon = 4;
	for (eslocal i = begin; i < end; i++) {
		for (size_t j = 0; j < elements[i]->coarseNodes(); j++) {
			metisNodes.push_back(localIndex(nodeProjection, elements[i]->node(j)));
		}
		if (nCommon > elements[i]->nCommon()) {
			nCommon = elements[i]->nCommon();
		}
	}
}

static std::vector<eslocal> continuousReorder(std::vector<Element*> &elements, size_t begin, size_t end)
{
	eslocal nCommon;
	std::vector<eslocal> nodeProjection, metisElements, metisNodes;

	METISMeshRepresentation(elements, begin, end, nodeProjection, metisElements, metisNodes, nCommon);

	eslocal ne = end - begin, nn = nodeProjection.size(), numflag = 0, *xAdj, *adjncy;
	checkMETISResult(METIS_MeshToDual(&ne, &nn, metisElements.data(), metisNodes.data(), &nCommon, &numflag, &xAdj, &adjncy));

	std::vector<eslocal> partPtrs;
	std::vector<Element*> eCopy(elements.begin() + begin, elements.begin() + end);

	partPtrs.push_back(begin);
	size_t back = begin;
	for (size_t e = 0; e < eCopy.size(); e++) {
		if (eCopy[e] != NULL) {
			std::vector<eslocal> stack;
			elements[back++] = eCopy[e];
			eCopy[e] = NULL;
			stack.push_back(e);
			while (stack.size()) {
				eslocal current = stack.back();
				stack.pop_back();
				for (size_t i = xAdj[current]; i < xAdj[current + 1]; i++) {
					if (eCopy[adjncy[i]] != NULL) {
						elements[back++] = eCopy[adjncy[i]];
						eCopy[adjncy[i]] = NULL;
						stack.push_back(adjncy[i]);
					}
				}
			}
			partPtrs.push_back(back);
		}
	}

	METIS_Free(xAdj);
	METIS_Free(adjncy);

	return partPtrs;
}

std::vector<eslocal> Mesh::getPartition(const std::vector<Element*> &elements, size_t begin, size_t end, eslocal parts) const
{
	if (parts == 1) {
		return std::vector<eslocal> (end - begin, 0);
	}
	eslocal nCommon;
	std::vector<eslocal> nodeProjection, metisElements, metisNodes;

	METISMeshRepresentation(elements, begin, end, nodeProjection, metisElements, metisNodes, nCommon);

	return computePartition(metisElements, metisNodes, end - begin, nodeProjection.size(), nCommon, parts);
}

std::vector<eslocal> Mesh::getPartition(size_t begin, size_t end, eslocal parts) const
{
	if (parts == 1) {
		return std::vector<eslocal> (end - begin, 0);
	}
	eslocal ncommon;
	std::vector<eslocal> e, n;

	// create array storing pointers to elements' nodes
	e.reserve(end - begin + 1);
	e.push_back(0);
	for (eslocal i = begin; i < end; i++) {
		e.push_back(e.back() + _elements[i]->coarseNodes());
	}

	// create array of nodes
	n.reserve(e.back());
	// number of common nodes to be neighbor
	ncommon = 4;
	for (eslocal i = begin; i < end; i++) {
		n.insert(n.end(), _elements[i]->indices(), _elements[i]->indices() + _elements[i]->coarseNodes());
		if (ncommon > _elements[i]->nCommon()) {
			ncommon = _elements[i]->nCommon();
		}
	}

	return computePartition(e, n, end - begin, _coordinates.clusterSize(), ncommon, parts);
}

static eslocal computeCenter(std::vector<std::vector<eslocal> > &neighbours)
{
	std::vector<eslocal> ia, ja;

	ia.reserve(neighbours.size() + 1);
	ia.push_back(0);
	for (size_t i = 0; i < neighbours.size(); i++) {
		std::sort(neighbours[i].begin(), neighbours[i].end());
		Esutils::removeDuplicity(neighbours[i]);
		ia.push_back(ia.back() + neighbours[i].size());
		ja.insert(ja.end(), neighbours[i].begin(), neighbours[i].end());
	}

	std::vector<float> a(ja.size(), 1);

	eslocal nSize = neighbours.size();
	std::vector<float> x(nSize, 1. / nSize), y(nSize);

	// Initial vector
	float last_l = nSize, l = 1;
	eslocal incr = 1;

	while (fabs((l - last_l) / l) > 1e-6) {
		mkl_cspblas_scsrsymv("U", &nSize, a.data(), ia.data(), ja.data(), x.data(), y.data());
		last_l = l;
		l = snrm2(&nSize, y.data(), &incr);
		cblas_sscal(nSize, 1 / l, y.data(), incr);
		x.swap(y);
	}

	return cblas_isamax(nSize, x.data(), incr);
}

eslocal Mesh::getCentralNode(const std::vector<Element*> &elements, size_t begin, size_t end, const std::vector<eslocal> &ePartition, eslocal subpart) const
{
	std::vector<eslocal> nMap;
	for (size_t e = begin; e < end; e++) {
		if (ePartition[e - begin] == subpart) {
			for (size_t n = 0; n < elements[e]->nodes(); n++) {
				nMap.push_back(elements[e]->node(n));
			}
		}
	}

	if (!nMap.size()) {
		return -1;
	}

	std::sort(nMap.begin(), nMap.end());
	Esutils::removeDuplicity(nMap);

	auto localIndex = [&] (eslocal index) {
		return std::lower_bound(nMap.begin(), nMap.end(), index) - nMap.begin();
	};

	std::vector<std::vector<eslocal> > neighbours(nMap.size());
	for (eslocal e = begin; e < end; e++) {
		if (ePartition[e - begin] == subpart) {
			for (size_t n = 0; n < elements[e]->nodes(); n++) {
				std::vector<eslocal> neigh = elements[e]->getNeighbours(n);
				for (size_t i = 0; i < neigh.size(); i++) {
					if (elements[e]->node(n) < neigh[i]) {
						neighbours[localIndex(elements[e]->node(n))].push_back(localIndex(neigh[i]));
					}
				}
			}
		}
	}

	return nMap[computeCenter(neighbours)];
}

eslocal Mesh::getCentralNode(eslocal begin, eslocal end, const std::vector<eslocal> &ePartition, eslocal part, eslocal subpart) const
{
	// Compute CSR format of symmetric adjacency matrix
	////////////////////////////////////////////////////////////////////////////
	std::vector<std::vector<eslocal> > neighbours(_coordinates.localSize(part));
	for (eslocal i = begin; i < end; i++) {
		if (ePartition[i - begin] == subpart) {
			for (size_t j = 0; j < _elements[i]->nodes(); j++) {
				std::vector<eslocal> neigh = _elements[i]->getNeighbours(j);
				for (size_t k = 0; k < neigh.size(); k++) {
					if (_elements[i]->node(j) < neigh[k]) {
						neighbours[_coordinates.localIndex(_elements[i]->node(j), part)].push_back(_coordinates.localIndex(neigh[k], part));
					}
				}
			}
		}
	}

	return _coordinates.clusterIndex(computeCenter(neighbours), part);
}

Mesh::~Mesh()
{
	for (size_t i = 0; i < _elements.size(); i++) {
		delete _elements[i];
	}

	for (size_t i = 0; i < _faces.size(); i++) {
		delete _faces[i];
	}

	for (size_t i = 0; i < _edges.size(); i++) {
		delete _edges[i];
	}

	for (size_t i = 0; i < _nodes.size(); i++) {
		delete _nodes[i];
	}

	for (size_t i = 0; i < _evaluators.size(); i++) {
		delete _evaluators[i];
	}
}

static bool isOuterFace(
		std::vector<std::vector<eslocal> > &nodesElements,
		std::vector<eslocal> &face)
{
	std::vector<eslocal> result(nodesElements[face[0]]);
	std::vector<eslocal>::iterator it = result.end();

	for (size_t i = 1; i < face.size(); i++) {
		std::vector<eslocal> tmp(result.begin(), it);
		it = std::set_intersection(tmp.begin(), tmp.end(),
				nodesElements[face[i]].begin(), nodesElements[face[i]].end(),
				result.begin());
		if (it - result.begin() == 1) {
			return true;
		}
	}

	return false;
}

void Mesh::getSurface(Mesh &surface) const
{
	// vector of faces in all parts
	std::vector<std::vector<std::vector<eslocal> > > faces(parts());
	// number of elements in all parts
	std::vector<size_t> elementsCount(parts(), 0);

	if (parts() < 1) {
		ESINFO(ERROR) << "Internal error: _partPtrs.size().";
	}

	cilk_for (size_t i = 0; i < parts(); i++) {
		// Compute nodes' adjacent elements
		std::vector<std::vector<eslocal> > nodesElements(_coordinates.localSize(i));
		for (eslocal j = _partPtrs[i]; j < _partPtrs[i + 1]; j++) {
			for (size_t k = 0; k < _elements[j]->nodes(); k++) {
				nodesElements[_elements[j]->node(k)].push_back(j);
			}
		}

		// compute number of elements and fill used nodes
		for (eslocal j = _partPtrs[i]; j < _partPtrs[i + 1]; j++) {
			for (size_t k = 0; k < _elements[j]->faces(); k++) {
				std::vector<eslocal> face(_elements[j]->face(k)->indices(), _elements[j]->face(k)->indices() + _elements[j]->face(k)->nodes());
				if (isOuterFace(nodesElements, face)) {
					for (size_t f = 0; f < face.size(); f++) {
						face[f] = _coordinates.clusterIndex(face[f], i);
					}
					faces[i].push_back(face);
					if (face.size() == 3) {
						elementsCount[i] += 1;
					}
					if (face.size() == 4) {
						elementsCount[i] += 2;
					}
				}
			}
		}
	}

	surface._coordinates = _coordinates;

	size_t count = 0;
	for (size_t i = 0; i + 1 < _partPtrs.size(); i++) {
		count += elementsCount[i];
	}

	surface._elements.reserve(count);
	surface._partPtrs.clear();
	surface._partPtrs.reserve(_partPtrs.size());
	eslocal params[6] = {0, 0, 0, 0, 0, 0};

	// create surface mesh
	surface._partPtrs.push_back(0); //(surface._elements.size());
	for (size_t i = 0; i + 1 < _partPtrs.size(); i++) {
		for (size_t j = 0; j < faces[i].size(); j++) {
			std::vector<eslocal> &face = faces[i][j];
			if (face.size() == 3) {
				surface._elements.push_back(new Triangle3(&face[0], params));
			}
			// divide square to triangles
			if (face.size() == 4) {
				size_t min = 0;
				for (size_t p = 1; p < 4; p++) {
					if (_coordinates[face[p]] < _coordinates[face[min]]) {
						min = p;
					}
				}
				if (min % 2 == 0) {
					surface._elements.push_back(new Triangle3(&face[0], params));
					face[1] = face[0];
					surface._elements.push_back(new Triangle3(&face[1], params));
				} else {
					surface._elements.push_back(new Triangle3(&face[1], params));
					face[2] = face[3];
					surface._elements.push_back(new Triangle3(&face[0], params));
				}
			}
		}
		surface._partPtrs.push_back(surface._elements.size());
	}

	surface.mapCoordinatesToDomains();
	surface._neighbours = _neighbours;
	surface._materials = _materials;
	for (size_t i = 0; i < _evaluators.size(); i++) {
		surface._evaluators.push_back(_evaluators[i]->copy());
	}
}

void Mesh::computeVolumeCorners(size_t number, bool vertices, bool edges, bool faces)
{
	if (parts() == 1 || (!vertices && !edges && !faces)) {
		return;
	}

	if (_corners.size()) {
		return;
	}

	computeFacesSharedByDomains();
	computeEdgesOnBordersOfFacesSharedByDomains();

	computeCornersOnEdges(number);
	if (faces) {
		computeCornersOnFaces(number);
	}
}

void Mesh::computePlaneCorners(size_t number, bool vertices, bool edges)
{
	if (parts() == 1 || (!vertices && !edges)) {
		return;
	}

	if (_corners.size()) {
		return;
	}

	computeEdgesSharedByDomains();
	computeCornersOnEdges(number);
}

template<typename MergeFunction>
static void uniqueWithMerge(std::vector<Element*> &elements, MergeFunction merge)
{
	if (elements.size() == 0) {
		return;
	}

	auto it = elements.begin();
	auto last = it;
	while (++it != elements.end()) {
		if (!(**last == **it)) {
			*(++last) = *it;
		} else { // Merge two edges
			merge(*last, *it);
			delete *it;
		}
	}
	elements.resize(++last - elements.begin());
}

template<typename MergeFunction>
static std::vector<Element*> mergeElements(size_t threads, std::vector<size_t> &distribution, std::vector<std::vector<Element*> > &elements, MergeFunction merge)
{
	std::vector<Element*> result;

	#pragma cilk grainsize = 1
	cilk_for (size_t t = 0; t < threads; t++) {
		std::sort(elements[t].begin(), elements[t].end(), [] (const Element* e1, const Element* e2) { return *e1 < *e2; });
		uniqueWithMerge(elements[t], merge);
	}

	std::vector<std::vector<Element*> > divided;
	std::vector<std::vector<Element*> > merged;

	divided.swap(elements);
	while (divided.size() > 1) {
		divided.resize(divided.size() + divided.size() % 2); // keep the size even
		merged.resize(divided.size() / 2);

		#pragma cilk grainsize = 1
		cilk_for (size_t t = 0; t < merged.size(); t++) {
			merged[t].resize(divided[2 * t].size() + divided[2 * t + 1].size());
			std::merge(
					divided[2 * t    ].begin(), divided[2 * t    ].end(),
					divided[2 * t + 1].begin(), divided[2 * t + 1].end(),
					merged[t].begin(), [] (const Element* e1, const Element* e2) { return *e1 < *e2; });
			uniqueWithMerge(merged[t], merge);
		}
		divided.swap(merged);
		merged.clear();
	}

	result.swap(divided[0]);
	return result;
}

void Mesh::fillEdgesFromElements(std::function<bool(const std::vector<Element*> &nodes, const Element* edge)> filter)
{
	size_t threads = config::env::CILK_NWORKERS;
	std::vector<size_t> distribution = Esutils::getDistribution(threads, _elements.size());

	std::vector<std::vector<Element*> > edges(threads);

	#pragma cilk grainsize = 1
	cilk_for (size_t t = 0; t < threads; t++) {
		for (size_t e = distribution[t]; e < distribution[t + 1]; e++) {
			_elements[e]->fillEdges();

			for (size_t i = 0; i < _elements[e]->edges(); i++) {
				Element* edge = _elements[e]->edge(i);
				if (edge->settings().size() || filter(_nodes, edge)) {
					edges[t].push_back(edge);
				} else {
					_elements[e]->setEdge(i, NULL);
					delete edge;
				}
			}
		}
	}

	_edges = mergeElements(threads, distribution, edges, [] (Element* e1, Element *e2) {
		e1->parentElements().insert(e1->parentElements().end(), e2->parentElements().begin(), e2->parentElements().end());
		for (size_t e = 0; e < e2->parentElements().size(); e++) {
			for (size_t i = 0; i < e2->parentElements()[e]->edges(); i++) {
				if (e2->parentElements()[0]->edge(i) != NULL && *(e1) == *(e2->parentElements()[e]->edge(i))) {
					e2->parentElements()[e]->setEdge(i, e1);
					break;
				}
			}
		}
	});

	mapEdgesToClusters();
	mapEdgesToDomains();
	fillParentEdgesToNodes();
}

void Mesh::fillFacesFromElements(std::function<bool(const std::vector<Element*> &nodes, const Element* face)> filter)
{
	size_t threads = config::env::CILK_NWORKERS;
	std::vector<size_t> distribution = Esutils::getDistribution(threads, _elements.size());

	std::vector<std::vector<Element*> > faces(threads);

	#pragma cilk grainsize = 1
	cilk_for (size_t t = 0; t < threads; t++) {
		for (size_t e = distribution[t]; e < distribution[t + 1]; e++) {
			_elements[e]->fillFaces();

			for (size_t f = 0; f < _elements[e]->faces(); f++) {
				Element* face = _elements[e]->face(f);
				if (face->settings().size() || filter(_nodes, face)) {
					faces[t].push_back(face);
				} else {
					_elements[e]->setFace(f, NULL);
					delete face;
				}
			}
		}
	}

	_faces = mergeElements(threads, distribution, faces, [] (Element* e1, Element *e2) {
		e1->parentElements().push_back(e2->parentElements().back()); // Face can be only between two elements
		for (size_t i = 0; i < e2->parentElements()[0]->faces(); i++) {
			if (e2->parentElements()[0]->face(i) != NULL && *e1 == *e2->parentElements()[0]->face(i)) {
				e2->parentElements()[0]->setFace(i, e1);
				break;
			}
		}
	});

	mapFacesToClusters();
	mapFacesToDomains();
	fillParentFacesToNodes();
}

void Mesh::fillNodesFromCoordinates()
{
	_nodes.reserve(_coordinates.clusterSize());
	for (size_t i = 0; i < _coordinates.clusterSize(); i++) {
		_nodes.push_back(new Node(i));
	}
}

void Mesh::fillParentElementsToNodes()
{
	cilk_for (size_t i = 0; i < _nodes.size(); i++) {
		_nodes[i]->parentElements().clear();
	}

	for (size_t e = 0; e < _elements.size(); e++) {
		for (size_t n = 0; n < _elements[e]->nodes(); n++) {
			_nodes[_elements[e]->node(n)]->parentElements().push_back(_elements[e]);
		}
	}

	cilk_for (size_t i = 0; i < _nodes.size(); i++) {
		std::sort(_nodes[i]->parentElements().begin(), _nodes[i]->parentElements().end());
	}
}

void APIMesh::fillParentElementsToDOFs(const std::vector<std::vector<eslocal> > &eDOFs)
{
	ESTEST(MANDATORY) << "Invalid number of recognized elements in API." << (eDOFs.size() != _elements.size() ? TEST_FAILED : TEST_PASSED);

	cilk_for (size_t i = 0; i < _DOFs.size(); i++) {
		_DOFs[i]->parentElements().clear();
	}

	for (size_t e = 0; e < eDOFs.size(); e++) {
		for (size_t d = 0; d < eDOFs[e].size(); d++) {
			_DOFs[eDOFs[e][d]]->parentElements().push_back(_elements[e]);
		}
	}

	cilk_for (size_t i = 0; i < _DOFs.size(); i++) {
		std::sort(_DOFs[i]->parentElements().begin(), _DOFs[i]->parentElements().end());
	}
}

void Mesh::fillParentFacesToNodes()
{
	cilk_for (size_t i = 0; i < _nodes.size(); i++) {
		_nodes[i]->parentFaces().clear();
	}

	for (size_t f = 0; f < _faces.size(); f++) {
		for (size_t n = 0; n < _faces[f]->nodes(); n++) {
			_nodes[_faces[f]->node(n)]->parentFaces().push_back(_faces[f]);
		}
	}

	cilk_for (size_t i = 0; i < _nodes.size(); i++) {
		std::sort(_nodes[i]->parentFaces().begin(), _nodes[i]->parentFaces().end());
	}
}

void Mesh::fillParentEdgesToNodes()
{
	cilk_for (size_t i = 0; i < _nodes.size(); i++) {
		_nodes[i]->parentEdges().clear();
	}

	for (size_t e = 0; e < _edges.size(); e++) {
		for (size_t n = 0; n < _edges[e]->nodes(); n++) {
			_nodes[_edges[e]->node(n)]->parentEdges().push_back(_edges[e]);
		}
	}

	cilk_for (size_t i = 0; i < _nodes.size(); i++) {
		std::sort(_nodes[i]->parentEdges().begin(), _nodes[i]->parentEdges().end());
	}
}

void Mesh::fillEdgesFromFaces(std::function<bool(const std::vector<Element*> &faces, const Element* edge)> filter)
{
	size_t threads = config::env::CILK_NWORKERS;
	std::vector<size_t> distribution = Esutils::getDistribution(threads, _faces.size());

	std::vector<std::vector<Element*> > edges(threads);

	#pragma cilk grainsize = 1
	cilk_for (size_t t = 0; t < threads; t++) {
		for (size_t e = distribution[t]; e < distribution[t + 1]; e++) {
			_faces[e]->fillEdges();

			for (size_t i = 0; i < _faces[e]->edges(); i++) {
				Element* edge = _faces[e]->edge(i);
				if (edge->settings().size() || filter(_nodes, edge)) {
					edges[t].push_back(edge);
				} else {
					_faces[e]->setEdge(i, NULL);
					delete edge;
				}
			}
		}
	}

	_edges = mergeElements(threads, distribution, edges, [] (Element* e1, Element *e2) {
		e1->parentElements().insert(e1->parentElements().end(), e2->parentElements().begin(), e2->parentElements().end());
		for (size_t e = 0; e < e2->parentElements().size(); e++) {
			for (size_t i = 0; i < e2->parentElements()[e]->edges(); i++) {
				if (e2->parentElements()[e]->edge(i) != NULL && *(e1) == *(e2->parentElements()[e]->edge(i))) {
					e2->parentElements()[e]->setEdge(i, e1);
					break;
				}
			}
		}
	});

	mapEdgesToClusters();
	mapEdgesToDomains();
	fillParentEdgesToNodes();
}

static Element* parentElement(const std::vector<Element*> &nodes, const Element *e)
{
	std::vector<Element*> intersection(nodes[e->node(e->nodes() - 1)]->parentElements()); // it is better to start from end (from mid points)
	auto it = intersection.end();

	for (size_t n = e->nodes() - 2; it - intersection.begin() > 1 &&  n < e->nodes(); n--) {
		std::vector<Element*> tmp(intersection.begin(), it);
		it = std::set_intersection(tmp.begin(), tmp.end(),
				nodes[e->node(n)]->parentElements().begin(), nodes[e->node(n)]->parentElements().end(),
				intersection.begin());
	}

	intersection.resize(it - intersection.begin());
	return intersection[0];
}

void Mesh::fillEdgesParents()
{
	size_t threads = config::env::CILK_NWORKERS;
	std::vector<size_t> distribution = Esutils::getDistribution(threads, _edges.size());

	#pragma cilk grainsize = 1
	cilk_for (size_t t = 0; t < threads; t++) {
		for (size_t e = distribution[t]; e < distribution[t + 1]; e++) {
			parentElement(_nodes, _edges[e])->setEdge(_edges[e]);
		}
	}
}

void Mesh::fillFacesParents()
{
	size_t threads = config::env::CILK_NWORKERS;
	std::vector<size_t> distribution = Esutils::getDistribution(threads, _faces.size());

	#pragma cilk grainsize = 1
	cilk_for (size_t t = 0; t < threads; t++) {
		for (size_t e = distribution[t]; e < distribution[t + 1]; e++) {
			parentElement(_nodes, _faces[e])->setEdge(_faces[e]);
		}
	}
}

void Mesh::computeFacesOfAllElements()
{
	fillFacesFromElements([] (const std::vector<Element*> &nodes, const Element* face) { return true; });
}


void Mesh::computeFacesOnDomainsSurface()
{
	fillFacesFromElements([] (const std::vector<Element*> &nodes, const Element *face) {
		std::vector<Element*> intersection(nodes[face->node(face->nodes() - 1)]->parentElements()); // it is better to start from end (from mid points)
		auto it = intersection.end();

		for (size_t n = face->nodes() - 2; it - intersection.begin() > 1 &&  n < face->nodes(); n--) {
			std::vector<Element*> tmp(intersection.begin(), it);
			it = std::set_intersection(tmp.begin(), tmp.end(),
					nodes[face->node(n)]->parentElements().begin(), nodes[face->node(n)]->parentElements().end(),
					intersection.begin());
		}

		intersection.resize(it - intersection.begin());
		if (intersection.size() == 1) {
			return true;
		}
		return intersection[0]->domains() != intersection[1]->domains();
	});
}

void Mesh::computeFacesSharedByDomains()
{
	fillFacesFromElements([] (const std::vector<Element*> &nodes, const Element *face) {
		std::vector<Element*> intersection(nodes[face->node(face->nodes() - 1)]->parentElements()); // it is better to start from end (from mid points)
		auto it = intersection.end();

		for (size_t n = face->nodes() - 2; it - intersection.begin() > 1 &&  n < face->nodes(); n--) {
			std::vector<Element*> tmp(intersection.begin(), it);
			it = std::set_intersection(tmp.begin(), tmp.end(),
					nodes[face->node(n)]->parentElements().begin(), nodes[face->node(n)]->parentElements().end(),
					intersection.begin());
		}

		intersection.resize(it - intersection.begin());
		if (intersection.size() == 1) {
			return false;
		}
		return intersection[0]->domains() != intersection[1]->domains();
	});
}

void Mesh::clearFacesWithoutSettings()
{
	size_t threads = config::env::CILK_NWORKERS;
	std::vector<size_t> distribution = Esutils::getDistribution(threads, _faces.size());

	#pragma cilk grainsize = 1
	cilk_for (size_t t = 0; t < threads; t++) {
		for (size_t f = distribution[t]; f < distribution[t + 1]; f++) {
			if (!_faces[f]->settings().size()) {
				for (size_t e = 0; e < _faces[f]->parentElements().size(); e++) {
					for (size_t i = 0; i < _faces[f]->parentElements()[e]->faces(); i++) {
						if (_faces[f]->parentElements()[e]->face(i) != NULL && *_faces[f] == *_faces[f]->parentElements()[e]->face(i)) {
							_faces[f]->parentElements()[e]->setFace(i, NULL);
						}
					}
				}
				delete _faces[f];
				_faces[f] = NULL;
			}
		}
	}

	size_t it = 0;
	for (size_t i = 0; i < _faces.size(); i++) {
		if (_faces[i] != NULL) {
			_faces[it++] = _faces[i];
		}
	}
	_faces.resize(it);
}

void Mesh::computeEdgesOfAllElements()
{
	fillEdgesFromElements([] (const std::vector<Element*> &nodes, const Element* edge) { return true; });
}

void Mesh::computeEdgesSharedByDomains()
{
	fillEdgesFromElements([] (const std::vector<Element*> &nodes, const Element *edge) {
		std::vector<Element*> intersection(nodes[edge->node(edge->nodes() - 1)]->parentElements()); // it is better to start from end (from mid points)
		auto it = intersection.end();

		for (size_t n = edge->nodes() - 2; it - intersection.begin() > 1 &&  n < edge->nodes(); n--) {
			std::vector<Element*> tmp(intersection.begin(), it);
			it = std::set_intersection(tmp.begin(), tmp.end(),
					nodes[edge->node(n)]->parentElements().begin(), nodes[edge->node(n)]->parentElements().end(),
					intersection.begin());
		}

		intersection.resize(it - intersection.begin());
		if (intersection.size() == 1) {
			return false;
		}
		return intersection[0]->domains() != intersection[1]->domains();
	});
}

void Mesh::computeEdgesOnBordersOfFacesSharedByDomains()
{
	fillEdgesFromFaces([] (const std::vector<Element*> &nodes, const Element* edge) {
		std::vector<Element*> intersection(nodes[edge->node(edge->nodes() - 1)]->parentFaces()); // it is better to start from end (from mid points)
		auto it = intersection.end();

		for (size_t n = edge->nodes() - 2; it - intersection.begin() > 1 &&  n < edge->nodes(); n--) {
			std::vector<Element*> tmp(intersection.begin(), it);
			it = std::set_intersection(tmp.begin(), tmp.end(),
					nodes[edge->node(n)]->parentFaces().begin(), nodes[edge->node(n)]->parentFaces().end(),
					intersection.begin());
		}

		intersection.resize(it - intersection.begin());
		if (intersection.size() == 1 && intersection[0]->domains().size() == 2) {
			return true;
		}

		for (size_t i = 1; i < intersection.size(); i++) {
			if (intersection[i]->domains() != intersection[i - 1]->domains()) {
				return true;
			}
		}
		return false;
	});
}

void Mesh::clearEdgesWithoutSettings()
{
	size_t threads = config::env::CILK_NWORKERS;
	std::vector<size_t> distribution = Esutils::getDistribution(threads, _edges.size());

	#pragma cilk grainsize = 1
	cilk_for (size_t t = 0; t < threads; t++) {
		for (size_t f = distribution[t]; f < distribution[t + 1]; f++) {
			if (!_edges[f]->settings().size()) {
				for (size_t e = 0; e < _edges[f]->parentElements().size(); e++) {
					for (size_t i = 0; i < _edges[f]->parentElements()[e]->edges(); i++) {
						if (_edges[f]->parentElements()[e]->edge(i) != NULL && *_edges[f] == *_edges[f]->parentElements()[e]->edge(i)) {
							_edges[f]->parentElements()[e]->setEdge(i, NULL);
						}
					}
				}
				delete _edges[f];
				_edges[f] = NULL;
			}
		}
	}

	size_t it = 0;
	for (size_t i = 0; i < _edges.size(); i++) {
		if (_edges[i] != NULL) {
			_edges[it++] = _edges[i];
		}
	}
	_edges.resize(it);
}

void Mesh::computeCornersOnEdges(size_t number)
{
	if (_edges.size() == 0) {
		ESINFO(ERROR) << "There are no edges for computation of corners.";
	}

	std::vector<Element*> edges;
	if (config::mesh::EDGE_CORNERS) {
		edges.reserve(_edges.size());
	}
	for (size_t e = 0; e < _edges.size(); e++) {
		if (_edges[e]->domains().size()) {
			if (config::mesh::VERTEX_CORNERS) {
				if (_nodes[_edges[e]->node(0)]->domains().size() > _nodes[_edges[e]->node(_edges[e]->coarseNodes() - 1)]->domains().size()) {
					_corners.push_back(_nodes[_edges[e]->node(0)]);
				}
				if (_nodes[_edges[e]->node(0)]->domains().size() < _nodes[_edges[e]->node(_edges[e]->coarseNodes() - 1)]->domains().size()) {
					_corners.push_back(_nodes[_edges[e]->node(_edges[e]->coarseNodes() - 1)]);
				}
			}
			if (config::mesh::EDGE_CORNERS) {
				edges.push_back(_edges[e]);
			}
		}
	}

	if (!edges.size()) {
		return;
	}

	std::sort(edges.begin(), edges.end(), [] (Element* e1, Element* e2) { return e1->domains() < e2->domains(); });

	std::vector<size_t> partPtrs;
	partPtrs.push_back(0);
	for (size_t i = 1; i < edges.size(); i++) {
		if (edges[i]->domains() != edges[i - 1]->domains()) {
			partPtrs.push_back(i);
		}
	}
	partPtrs.push_back(edges.size());

	std::vector<std::vector<Element*> > corners(partPtrs.size() - 1);
	cilk_for (size_t p = 0; p < partPtrs.size() - 1; p++) {
		std::vector<eslocal> subPartPtrs = continuousReorder(edges, partPtrs[p], partPtrs[p + 1]);
		for (size_t sp = 0; sp < subPartPtrs.size() - 1; sp++) {
			std::vector<eslocal> nodes;
			for (size_t e = subPartPtrs[sp]; e < subPartPtrs[sp + 1]; e++) {
				nodes.push_back(edges[e]->node(0));
				nodes.push_back(edges[e]->node(edges[e]->coarseNodes() - 1));
			}
			std::sort(nodes.begin(), nodes.end());
			bool cycle = nodes.size() % 2 == 0;
			for (size_t i = 0; cycle && i < nodes.size(); i += 2) {
				cycle = nodes[i + 1] == nodes[i];
			}
			if (cycle && nodes.size() < 12) {
				Esutils::removeDuplicity(nodes);
				for (size_t i = 0; i < nodes.size() && i < 4; i++) {
					corners[p].push_back(_nodes[nodes[i]]);
				}
				continue;
			}
			if (!cycle && !number) {
				continue;
			}
			std::vector<eslocal> ePartition = getPartition(edges, subPartPtrs[sp], subPartPtrs[sp + 1], cycle ? 4 : number);
			for (size_t c = 0; c < (cycle ? 4 : number); c++) {
				eslocal center = getCentralNode(edges, subPartPtrs[sp], subPartPtrs[sp + 1], ePartition, c);
				if (center > -1) {
					corners[p].push_back(_nodes[center]);
				}
			}
		}
	}

	for (size_t p = 0; p < corners.size(); p++) {
		_corners.insert(_corners.end(), corners[p].begin(), corners[p].end());
	}
	std::sort(_corners.begin(), _corners.end());
	Esutils::removeDuplicity(_corners);
}

void Mesh::computeCornersOnFaces(size_t number)
{
	ESINFO(GLOBAL_ERROR) << "Corners in faces are not implemented.";
}

static void setCluster(Element* &element, std::vector<Element*> &nodes)
{
	std::vector<eslocal> intersection = nodes[element->node(0)]->clusters();
	for (size_t i = 1; i < element->nodes(); i++) {
		auto tmp(intersection);
		auto it = std::set_intersection(
				nodes[element->node(i)]->clusters().begin(), nodes[element->node(i)]->clusters().end(),
				tmp.begin(), tmp.end(), intersection.begin());
		intersection.resize(it - intersection.begin());
	}
	element->clusters() = intersection;
}

void Mesh::mapFacesToClusters()
{
	size_t threads = config::env::CILK_NWORKERS;
	std::vector<size_t> distribution = Esutils::getDistribution(threads, _faces.size());

	#pragma cilk grainsize = 1
	cilk_for (size_t t = 0; t < threads; t++) {
		for (size_t f = distribution[t]; f < distribution[t + 1]; f++) {
			if (_faces[f]->parentElements().size() == 1) { // Only faces with one element can have more clusters
				if (std::all_of(_faces[f]->indices(), _faces[f]->indices() + _faces[f]->coarseNodes(), [&] (eslocal i) { return _nodes[i]->clusters().size() > 1; })) {
					setCluster(_faces[f], _nodes);
				}
			}
		}
	}
}

void Mesh::mapEdgesToClusters()
{
	size_t threads = config::env::CILK_NWORKERS;
	std::vector<size_t> distribution = Esutils::getDistribution(threads, _edges.size());

	#pragma cilk grainsize = 1
	cilk_for (size_t t = 0; t < threads; t++) {
		for (size_t e = distribution[t]; e < distribution[t + 1]; e++) {
			if (std::all_of(_edges[e]->indices(), _edges[e]->indices() + _edges[e]->coarseNodes(), [&] (eslocal i) { return _nodes[i]->clusters().size() > 1; })) {
				setCluster(_edges[e], _nodes);
			}
		}
	}
}

static void assignDomains(std::vector<Element*> &elements)
{
	size_t threads = config::env::CILK_NWORKERS;
	std::vector<size_t> distribution = Esutils::getDistribution(threads, elements.size());

	#pragma cilk grainsize = 1
	cilk_for (size_t t = 0; t < threads; t++) {
		for (size_t i = distribution[t]; i < distribution[t + 1]; i++) {
			elements[i]->domains().clear();
			for (size_t e = 0; e < elements[i]->parentElements().size(); e++) {
				elements[i]->domains().insert(elements[i]->domains().end(), elements[i]->parentElements()[e]->domains().begin(), elements[i]->parentElements()[e]->domains().end());
			}
			std::sort(elements[i]->domains().begin(), elements[i]->domains().end());
			Esutils::removeDuplicity(elements[i]->domains());
		}
	}
}

void Mesh::mapElementsToDomains()
{
	cilk_for (size_t p = 0; p < parts(); p++) {
		for (size_t e = _partPtrs[p]; e < _partPtrs[p + 1]; e++) {
			_elements[e]->domains().clear();
			_elements[e]->domains().push_back(p);
		}
	}
}

void Mesh::mapFacesToDomains()
{
	assignDomains(_faces);
}

void Mesh::mapEdgesToDomains()
{
	assignDomains(_edges);
}

void Mesh::mapNodesToDomains()
{
	assignDomains(_nodes);
}

void APIMesh::mapDOFsToDomains()
{
	assignDomains(_DOFs);
}

static void setDOFsIndices(
		std::vector<Element*> &elements,
		size_t parts,
		const std::vector<Property> &DOFs,
		const std::vector<size_t> &offsets,
		const std::vector<std::vector<std::vector<size_t> > > &threadsOffsets)
{
	size_t threads = config::env::CILK_NWORKERS;
	std::vector<size_t> distribution = Esutils::getDistribution(threads, elements.size());

	#pragma cilk grainsize = 1
	cilk_for (size_t t = 0; t < threads; t++) {
		std::vector<std::vector<size_t> > counters(parts);
		for (size_t p = 0; p < parts; p++) {
			switch (config::assembler::DOFS_ORDER) {
			case config::assembler::DOFS_ORDERalternative::GROUP_ELEMENTS:
				counters[p].resize(1, 0);
				for (size_t dof = 0; dof < DOFs.size(); dof++) {
					counters[p][0] += threadsOffsets[p][dof][t];
				}
				counters[p][0] += offsets[p];
				break;
			case config::assembler::DOFS_ORDERalternative::GROUP_DOFS:
				counters[p].resize(3, offsets[p]);
				counters[p][0] += threadsOffsets[p][0][t];
				counters[p][1] += threadsOffsets[p][1][t] + threadsOffsets[p][0][threads];
				counters[p][2] += threadsOffsets[p][2][t] + threadsOffsets[p][0][threads] + threadsOffsets[p][1][threads];
				break;
			}
		}

		for (size_t i = distribution[t]; i < distribution[t + 1]; i++) {
			for (size_t d = 0; d < elements[i]->domains().size(); d++) {

				for (size_t dof = 0; dof < DOFs.size(); dof++) {
					if (elements[i]->DOFsIndices()[d * DOFs.size() + dof] == 1) {

						switch (config::assembler::DOFS_ORDER) {
						case config::assembler::DOFS_ORDERalternative::GROUP_ELEMENTS:
							elements[i]->DOFsIndices()[d * DOFs.size() + dof] = counters[elements[i]->domains()[d]][0]++;
							break;
						case config::assembler::DOFS_ORDERalternative::GROUP_DOFS:
							elements[i]->DOFsIndices()[d * DOFs.size() + dof] = counters[elements[i]->domains()[d]][dof]++;
							break;
						}

					}
				}
			}
		}
	}
}

std::vector<size_t> Mesh::assignVariousDOFsIndicesToNodes(const std::vector<size_t> &offsets, const std::vector<Property> &DOFs)
{
	auto findDOF = [&] (const std::vector<Property> &DOFs, size_t &added, std::vector<bool> &addDOF) {
		for (size_t dof = 0; dof < DOFs.size(); dof++) {
			if (!addDOF[dof] && std::find(DOFs.begin(), DOFs.end(), DOFs[dof]) != DOFs.end()) {
				added++;
				addDOF[dof] = true;
			}
		}
	};

	auto fillDOFs = [&] (Element *node, eslocal i, eslocal domain, std::vector<bool> &addDOF) {
		for (size_t e = 0, added = 0; e < node->parentElements().size() && added < DOFs.size(); e++) {
			const Element *el = node->parentElements()[e];
			if (el->domains()[0] == domain) {
				if (std::find(el->indices(), el->indices() + el->coarseNodes(), i) == el->indices() + el->coarseNodes()) {
					findDOF(el->midPointDOFs(), added, addDOF);
				} else {
					findDOF(el->pointDOFs(), added, addDOF);
				}
			}
		}
	};


	size_t threads = config::env::CILK_NWORKERS;
	std::vector<size_t> distribution = Esutils::getDistribution(threads, _nodes.size());

	// domains x DOFs x (threads + 1)
	std::vector<std::vector<std::vector<size_t> > > threadsOffsets(parts(), std::vector<std::vector<size_t> >(DOFs.size(), std::vector<size_t>(threads + 1)));

	#pragma cilk grainsize = 1
	cilk_for (size_t t = 0; t < threads; t++) {
		std::vector<std::vector<size_t> > threadOffset(parts(), std::vector<size_t>(DOFs.size(), 0));
		for (size_t i = distribution[t]; i < distribution[t + 1]; i++) {
			_nodes[i]->DOFsIndices().resize(DOFs.size() * _nodes[i]->domains().size(), -1);

			for (size_t d = 0; d < _nodes[i]->domains().size(); d++) {
				std::vector<bool> addDOF(DOFs.size(), false);

				fillDOFs(_nodes[i], i, _nodes[i]->domains()[d], addDOF);

				for (size_t dof = 0; dof < DOFs.size(); dof++) {
					if (addDOF[dof]) {
						_nodes[i]->DOFsIndices()[d * DOFs.size() + dof] = 1;
						threadOffset[_nodes[i]->domains()[d]][dof]++;
					}
				}
			}
		}

		for (size_t p = 0; p < parts(); p++) {
			for (size_t d = 0; d < DOFs.size(); d++) {
				threadsOffsets[p][d][t] = threadOffset[p][d];
			}
		}
	}

	std::vector<size_t> sizes(offsets);
	for (size_t p = 0; p < parts(); p++) {
		for (size_t d = 0; d < DOFs.size(); d++) {
			sizes[p] += Esutils::sizesToOffsets(threadsOffsets[p][d]);
		}
	}

	setDOFsIndices(_nodes, parts(), DOFs, offsets, threadsOffsets);

	return sizes;
}


static std::vector<size_t> fillUniformDOFs(
		std::vector<Element*> &elements,
		size_t parts,
		const std::vector<Property> &DOFs,
		const std::vector<size_t> &offsets)
{
	size_t threads = config::env::CILK_NWORKERS;
	std::vector<size_t> distribution = Esutils::getDistribution(threads, elements.size());

	// domains x DOF x (threads + 1)
	std::vector<std::vector<std::vector<size_t> > > threadsOffsets(parts, std::vector<std::vector<size_t> >(DOFs.size(), std::vector<size_t>(threads + 1)));

	#pragma cilk grainsize = 1
	cilk_for (size_t t = 0; t < threads; t++) {
		std::vector<size_t> threadOffset(parts, 0);
		for (size_t i = distribution[t]; i < distribution[t + 1]; i++) {
			elements[i]->DOFsIndices().resize(DOFs.size() * elements[i]->domains().size(), 1);

			for (size_t d = 0; d < elements[i]->domains().size(); d++) {
				threadOffset[elements[i]->domains()[d]]++;
			}
		}

		for (size_t p = 0; p < parts; p++) {
			for (size_t d = 0; d < DOFs.size(); d++) {
				threadsOffsets[p][d][t] = threadOffset[p];
			}
		}
	}

	std::vector<size_t> sizes(offsets);
	for (size_t p = 0; p < parts; p++) {
		for (size_t d = 0; d < DOFs.size(); d++) {
			sizes[p] += Esutils::sizesToOffsets(threadsOffsets[p][d]);
		}
	}

	setDOFsIndices(elements, parts, DOFs, offsets, threadsOffsets);

	return sizes;
}

std::vector<size_t> Mesh::assignUniformDOFsIndicesToNodes(const std::vector<size_t> &offsets, const std::vector<Property> &DOFs)
{
	return fillUniformDOFs(_nodes, parts(), DOFs, offsets);
}

std::vector<size_t> Mesh::assignUniformDOFsIndicesToEdges(const std::vector<size_t> &offsets, const std::vector<Property> &DOFs)
{
	return fillUniformDOFs(_edges, parts(), DOFs, offsets);
}

std::vector<size_t> Mesh::assignUniformDOFsIndicesToFaces(const std::vector<size_t> &offsets, const std::vector<Property> &DOFs)
{
	return fillUniformDOFs(_faces, parts(), DOFs, offsets);
}

std::vector<size_t> Mesh::assignUniformDOFsIndicesToElements(const std::vector<size_t> &offsets, const std::vector<Property> &DOFs)
{
	return fillUniformDOFs(_elements, parts(), DOFs, offsets);
}

std::vector<size_t> APIMesh::distributeDOFsToDomains(const std::vector<size_t> &offsets)
{
	return fillUniformDOFs(_DOFs, parts(), { Property::UNKNOWN }, offsets);
}

static void computeDOFsCounters(std::vector<Element*> &elements, const std::vector<Property> &DOFs, const Mesh &mesh)
{
	std::vector<int> neighbours = mesh.neighbours();
	neighbours.push_back(config::env::MPIrank);
	std::sort(neighbours.begin(), neighbours.end());

	auto n2i = [ & ] (size_t neighbour) {
		return std::lower_bound(neighbours.begin(), neighbours.end(), neighbour) - neighbours.begin();
	};

	size_t threads = config::env::CILK_NWORKERS;
	std::vector<size_t> distribution = Esutils::getDistribution(threads, elements.size());

	// threads x neighbour x data
	std::vector<std::vector<std::vector<esglobal> > > sBuffer(threads, std::vector<std::vector<esglobal> >(neighbours.size()));
	// neighbour x data
	std::vector<std::vector<esglobal> > rBuffer(neighbours.size());

	// Compute send buffers
	#pragma cilk grainsize = 1
	cilk_for (size_t t = 0; t < threads; t++) {
		for (size_t e = distribution[t]; e < distribution[t + 1]; e++) {

			elements[e]->DOFsDomainsCounters().resize(DOFs.size() * elements[e]->clusters().size(), -1);
			size_t cluster = std::lower_bound(elements[e]->clusters().begin(), elements[e]->clusters().end(), config::env::MPIrank) - elements[e]->clusters().begin();
			for (size_t i = 0; i < DOFs.size(); i++) {
				elements[e]->DOFsDomainsCounters()[cluster * DOFs.size() + i] = elements[e]->numberOfLocalDomainsWithDOF(i);
			}
			if (elements[e]->clusters().size() > 1) {
				for (auto c = elements[e]->clusters().begin(); c != elements[e]->clusters().end(); ++c) {
					if (*c == config::env::MPIrank) {
						continue;
					}

					sBuffer[t][n2i(*c)].push_back(elements[e]->vtkCode());
					for (size_t n = 0; n < elements[e]->coarseNodes(); n++) {
						sBuffer[t][n2i(*c)].push_back(mesh.coordinates().globalIndex(elements[e]->node(n)));
					}

					for (size_t i = 0; i < DOFs.size(); i++) {
						sBuffer[t][n2i(*c)].push_back(elements[e]->DOFsDomainsCounters()[cluster * DOFs.size() + i]);
					}
				}
			}

		}
	}

	for (size_t t = 1; t < threads; t++) {
		for (size_t n = 0; n < neighbours.size(); n++) {
			sBuffer[0][n].insert(sBuffer[0][n].end(), sBuffer[t][n].begin(), sBuffer[t][n].end());
		}
	}

	std::vector<MPI_Request> req(neighbours.size());
	for (size_t n = 0; n < neighbours.size(); n++) {
		MPI_Isend(sBuffer[0][n].data(), sizeof(esglobal) * sBuffer[0][n].size(), MPI_BYTE, neighbours[n], 0, MPI_COMM_WORLD, req.data() + n);
	}

	int flag, counter = 0;
	MPI_Status status;
	while (counter < neighbours.size()) {
		MPI_Iprobe(MPI_ANY_SOURCE, 0, MPI_COMM_WORLD, &flag, &status);
		if (flag) {
			int count;
			MPI_Get_count(&status, MPI_BYTE, &count);
			rBuffer[n2i(status.MPI_SOURCE)].resize(count / sizeof(esglobal));
			MPI_Recv(rBuffer[n2i(status.MPI_SOURCE)].data(), count, MPI_BYTE, status.MPI_SOURCE, 0, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
			counter++;
		}
	}

	MPI_Waitall(neighbours.size(), req.data(), MPI_STATUSES_IGNORE);
	// All data are exchanged

	// Create list of neighbours elements
	std::vector<std::vector<Element*> > nElements(neighbours.size());
	#pragma cilk grainsize = 1
	cilk_for (size_t n = 0; n < neighbours.size(); n++) {
		size_t p = 0;
		while (p + 1 < rBuffer[n].size()) {
			switch (rBuffer[n][p++]) {
			case NodeVTKCode:
				nElements[n].push_back(new Node(rBuffer[n][p]));
				break;
			case Square4VTKCode:
				nElements[n].push_back(new Square4(&rBuffer[n][p]));
				break;
			case Square8VTKCode:
				nElements[n].push_back(new Square8(&rBuffer[n][p]));
				break;
			case Triangle3VTKCode:
				nElements[n].push_back(new Triangle3(&rBuffer[n][p]));
				break;
			case Triangle6VTKCode:
				nElements[n].push_back(new Triangle6(&rBuffer[n][p]));
				break;
			default:
				// Volume elements are never exchanged
				ESINFO(GLOBAL_ERROR) << "Unknown neighbour element";
			}
			p += nElements[n].back()->nodes();
			for (size_t i = 0; i < nElements[n].back()->nodes(); i++) {
				nElements[n].back()->node(i) = mesh.coordinates().clusterIndex(nElements[n].back()->node(i));
			}

			nElements[n].back()->DOFsDomainsCounters() = std::vector<eslocal>(&rBuffer[n][p], &rBuffer[n][p] + DOFs.size());
			p += DOFs.size();
		}
	}

	// TODO: parallelization
	for (size_t n = 0; n < neighbours.size(); n++) {
		for (size_t e = 0; e < nElements[n].size(); e++) {
			auto it = std::lower_bound(elements.begin(), elements.end(), nElements[n][e], [&] (Element *el1, Element *el2) { return *el1 < *el2; });
			if (it != elements.end() && **it == *(nElements[n][e])) {
				size_t cluster = std::lower_bound((*it)->clusters().begin(), (*it)->clusters().end(), neighbours[n]) - (*it)->clusters().begin();
				for (size_t dof = 0; dof < DOFs.size(); dof++) {
					(*it)->DOFsDomainsCounters()[cluster * DOFs.size() + dof] = nElements[n][e]->DOFsDomainsCounters()[dof];
				}
			}
		}
	}

	#pragma cilk grainsize = 1
	cilk_for (size_t n = 0; n < neighbours.size(); n++) {
		for (size_t e = 0; e < nElements[n].size(); e++) {
			delete nElements[n][e];
		}
	}

	// Remove elements that are not in both clusters
	#pragma cilk grainsize = 1
	cilk_for (size_t t = 0; t < threads; t++) {
		for (size_t e = distribution[t]; e < distribution[t + 1]; e++) {

			if (elements[e]->clusters().size() > 1) {
				std::vector<eslocal> &counters = elements[e]->DOFsDomainsCounters();
				for (size_t c = 0; c < elements[e]->clusters().size(); c++) {
					if (std::all_of(counters.begin() + c * DOFs.size(), counters.begin() + (c + 1) * DOFs.size(), [] (eslocal &c) { return c == -1; })) {
						counters.erase(counters.begin() + c * DOFs.size(), counters.begin() + (c + 1) * DOFs.size());
						elements[e]->clusters().erase(elements[e]->clusters().begin() + c--);
					}
				}
			}

		}
	}
}


void Mesh::computeNodesDOFsCounters(const std::vector<Property> &DOFs)
{
	computeDOFsCounters(_nodes, DOFs, *this);
}

void Mesh::computeEdgesDOFsCounters(const std::vector<Property> &DOFs)
{
	computeDOFsCounters(_edges, DOFs, *this);
}

void Mesh::computeFacesDOFsCounters(const std::vector<Property> &DOFs)
{
	computeDOFsCounters(_faces, DOFs, *this);
}

void APIMesh::computeDOFsDOFsCounters()
{
	computeDOFsCounters(_DOFs, { Property::UNKNOWN }, *this);
}

void Mesh::mapCoordinatesToDomains()
{
	_coordinates._clusterIndex.clear();
	_coordinates._clusterIndex.resize(parts());

	for (size_t p = 0; p < parts(); p++) {
		std::vector<eslocal> l2g;
		for (eslocal e = _partPtrs[p]; e < _partPtrs[p + 1]; e++) {
			l2g.insert(l2g.end(), _elements[e]->indices(), _elements[e]->indices() + _elements[e]->nodes());
		}

		std::sort(l2g.begin(), l2g.end());
		Esutils::removeDuplicity(l2g);

		_coordinates._clusterIndex[p] = l2g;
	}
}

}


