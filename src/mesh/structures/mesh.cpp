
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
		_partPtrs.resize(parts + 1);
		_partPtrs[0] = 0;
		_partPtrs[1] = _elements.size();
		mapElementsToDomains();
	} else {
		_partPtrs.resize(parts + 1);
		_partPtrs[0] = 0;

		eslocal *ePartition = getPartition(0, _elements.size(), parts);

		_partPtrs = std::vector<eslocal>(parts + 1, 0);
		for (size_t i = 0; i < _elements.size(); i++) {
			_elements[i]->domains().clear();
			_elements[i]->domains().push_back(ePartition[i]);
			_partPtrs[ePartition[i]]++;
		}

		std::sort(_elements.begin(), _elements.end(), [] (const Element* e1, const Element* e2) { return e1->domains()[0] < e2->domains()[0]; });
		ESTEST(MANDATORY) << "subdomain without element" << (std::any_of(_partPtrs.begin(), _partPtrs.end() - 1, [] (eslocal size) { return size == 0; }) ? TEST_FAILED : TEST_PASSED);
		Esutils::sizesToOffsets(_partPtrs);

		delete[] ePartition;
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
	ESINFO(GLOBAL_ERROR) << "Fix API partition function";
	if (parts == 1 && this->parts() == 1) {
		_partPtrs.resize(parts + 1);
		_partPtrs[0] = 0;
		_partPtrs[1] = _elements.size();
		mapElementsToDomains();
		return;
	}

	eslocal *ePartition = getPartition(0, _elements.size(), parts);

	_partPtrs = std::vector<eslocal>(parts + 1, 0);
	for (size_t i = 0; i < _elements.size(); i++) {
		_elements[i]->domains().push_back(ePartition[i]);
		_partPtrs[ePartition[i]]++;
	}

	std::sort(_elements.begin(), _elements.end(), [] (const Element* e1, const Element* e2) { return e1->domains()[0] < e2->domains()[0]; });
	ESTEST(MANDATORY) << "subdomain without element" << (std::any_of(_partPtrs.begin(), _partPtrs.end() - 1, [] (eslocal size) { return size == 0; }) ? TEST_FAILED : TEST_PASSED);
	Esutils::sizesToOffsets(_partPtrs);

	delete[] ePartition;
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
		eslocal *eSubPartition = getPartition(_partPtrs[part], _partPtrs[part + 1], points);

		for (eslocal j = 0; j < points; j++) {
			fixPoints[j] = getCentralNode(_partPtrs[part], _partPtrs[part + 1], eSubPartition, part, j);
		}
		std::sort(fixPoints.begin(), fixPoints.end());

		// Remove the same points
		Esutils::removeDuplicity(fixPoints);

		delete[] eSubPartition;

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

eslocal* Mesh::getPartition(eslocal first, eslocal last, eslocal parts) const
{
	if (parts == 1) {
		eslocal *ePartition = new eslocal[last - first];
		for (eslocal i = first; i < last; i++) {
			ePartition[i - first] = 0;
		}
		return ePartition;
	}
	// INPUTS
	eslocal ncommon, eSize, nSize, *e, *n, options[METIS_NOPTIONS];

	// OUTPUTS
	eslocal objval, *ePartition, *nPartition;

	// FILL INPUT VARIABLES
	////////////////////////////////////////////////////////////////////////////

	// There is probably BUG in METIS numbering or I do not understand documentation.
	// The solution is increase the size of 'nodesCount' and keep the default numbering
	//options[METIS_OPTION_NUMBERING] = coordinates.getNumbering();
	METIS_SetDefaultOptions(options);
	//TODO:
	options[METIS_OPTION_CONTIG]  = 1;
//	options[METIS_OPTION_MINCONN] = 1;
//	options[METIS_OPTION_NITER]   = 20;
//	options[METIS_PTYPE_KWAY]     = 1;

	eSize = last - first;
	nSize = _coordinates.clusterSize();

	// create array storing pointers to elements' nodes
	e = new eslocal[eSize + 1];
	e[0] = 0;
	for (eslocal i = first, index = 0; i < last; i++, index++) {
		e[index + 1] = e[index] + _elements[i]->coarseNodes();
	}

	// create array of nodes
	n = new eslocal[e[eSize]];
	// number of common nodes to be neighbor
	ncommon = 4;
	for (eslocal i = first, index = 0; i < last; i++, index++) {
		const Element* el = _elements[i];
		memcpy(n + e[index], el->indices(), el->coarseNodes() * sizeof(eslocal));
		if (ncommon > _elements[i]->nCommon()) {
			ncommon = _elements[i]->nCommon();
		}
	}

	// PREPARE OUTPUT VARIABLES
	////////////////////////////////////////////////////////////////////////////
	ePartition = new eslocal[eSize];
	nPartition = new eslocal[nSize];

	eslocal result = METIS_PartMeshDual(
			&eSize,
			&nSize,
			e,
			n,
			NULL,		// weights of nodes
			NULL,		// size of nodes
			&ncommon,
			&parts,
			NULL,		// weights of parts
			options,
			&objval,
			ePartition,
			nPartition);
	checkMETISResult(result);

	delete[] e;
	delete[] n;
	delete[] nPartition;

	return ePartition;
}

eslocal Mesh::getCentralNode(
		eslocal first,
		eslocal last,
		eslocal *ePartition,
		eslocal part,
		eslocal subpart) const
{
	// Compute CSR format of symmetric adjacency matrix
	////////////////////////////////////////////////////////////////////////////
	std::vector<std::set<eslocal> > neighbours(_coordinates.localSize(part));
	for (eslocal i = first, index = 0; i < last; i++, index++) {
		if (ePartition[index] == subpart) {
			for (size_t j = 0; j < _elements[i]->nodes(); j++) {
				std::vector<eslocal> neigh = _elements[i]->getNeighbours(j);
				for (size_t k = 0; k < neigh.size(); k++) {
					if (_elements[i]->node(j) < neigh[k]) {
						neighbours[_coordinates.localIndex(_elements[i]->node(j), part)].insert(_coordinates.localIndex(neigh[k], part));
					}
				}
			}
		}
	}
	eslocal nonZeroValues = 0;
	for (size_t i = 0; i < neighbours.size(); i++) {
		nonZeroValues += neighbours[i].size();
	}

	float *a;
	eslocal *ia, *ja;
	a = new float[nonZeroValues];
	ia = new eslocal[_coordinates.localSize(part) + 1];
	ja = new eslocal[nonZeroValues];

	eslocal i = 0;
	std::set<eslocal>::const_iterator it;
	for (size_t n = 0; n < neighbours.size(); n++) {
		std::set<eslocal> &values = neighbours[n];
		ia[n] = i;
		for (it = values.begin(); it != values.end(); ++it, i++) {
			a[i] = 1;
			ja[i] = *it;
		}
	}
	ia[neighbours.size()] = i;

	eslocal nSize = _coordinates.localSize(part);
	float *x, *y, *swap;
	x = new float[nSize];
	y = new float[nSize];

	// Initial vector
	for (eslocal xi = 0; xi < nSize; xi++) {
		x[xi] = 1. / nSize;
	}
	float last_l = nSize, l = 1;
	eslocal incr = 1;

	while (fabs((l - last_l) / l) > 1e-6) {
		mkl_cspblas_scsrsymv("U", &nSize, a, ia, ja, x, y);
		last_l = l;
		l = snrm2(&nSize, y, &incr);
		cblas_sscal(nSize, 1 / l, y, incr);
		swap = x;
		x = y;
		y = swap;
	}
	eslocal result = cblas_isamax(nSize, x, incr);

	delete[] a;
	delete[] ia;
	delete[] ja;
	delete[] x;
	delete[] y;

	return _coordinates.clusterIndex(result, part);
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

void Mesh::makePartContinuous(size_t part)
{
	size_t begin = _partPtrs[part];
	size_t end = _partPtrs[part + 1];

	eslocal ne, nn, *ePtr, *eInd, nCommon, numflag, *xAdj, *adjncy;

	ne = end - begin;
	nn = _coordinates.clusterSize();

	ePtr = new eslocal[ne + 1];
	ePtr[0] = 0;
	for (size_t e = begin, i = 1; e < end; e++, i++) {
		ePtr[i] = ePtr[i - 1] + _elements[e]->coarseNodes();
	}

	nCommon = 4;
	eInd = new eslocal[ePtr[ne]];
	for (size_t e = begin, i = 0; e < end; e++, i++) {
		const Element* el = _elements[e];
		memcpy(eInd + ePtr[i], el->indices(), el->coarseNodes() * sizeof(eslocal));
		if (nCommon > _elements[i]->nCommon()) {
			nCommon = _elements[i]->nCommon();
		}
	}

	numflag = 0;
	checkMETISResult(METIS_MeshToDual(&ne, &nn, ePtr, eInd, &nCommon, &numflag, &xAdj, &adjncy));

	std::vector<int> color(ne, 0);

	std::function<void(eslocal, int)> traverse = [&] (eslocal element, int eColor) {
		if (color[element] > 0) {
			return;
		}
		color[element] = eColor;
		for (eslocal e = xAdj[element]; e < xAdj[element + 1]; e++) {
			traverse(adjncy[e], eColor);
		}
	};

	int colors = 0;
	for (size_t i = 0; i < color.size(); i++) {
		if (color[i] == 0) {
			traverse(i, ++colors);
		}
	}

	delete[] ePtr;
	delete[] eInd;
	METIS_Free(xAdj);
	METIS_Free(adjncy);

	if (colors == 1) {
		return;
	}

	std::vector<Element*> tmp(_elements.begin() + begin, _elements.begin() + end);
	size_t p = begin;
	for (int c = 1; c <= colors; c++) {
		for (size_t e = 0; e < tmp.size(); e++) {
			if (color[e] == c) {
				_elements[p++] = tmp[e];
			}
		}
		if (c < colors) {
			_partPtrs.insert(_partPtrs.begin() + part + c, p);
		}
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

void Mesh::fillNodesFromElements()
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

	std::vector<Element*> edges = _edges;
	std::sort(edges.begin(), edges.end(), [] (Element* e1, Element* e2) { return e1->domains() < e2->domains(); });

	auto is_cycle = [] (std::vector<eslocal> &nodes) {
		std::sort(nodes.begin(), nodes.end());
		size_t fullSize = nodes.size();
		Esutils::removeDuplicity(nodes);
		return fullSize / 2 == nodes.size();
	};

	auto next = [&] (Element* &edge, Element* &node) {
		edge = edge == node->parentEdges()[0] ? node->parentEdges()[1] : node->parentEdges()[0];
		node = node == _nodes[edge->node(0)] ? _nodes[edge->node(edge->nodes() - 1)] : _nodes[edge->node(0)];

	};

	auto fixCycle = [&] (size_t begin, size_t end) {
		Element *edge = edges[begin], *node = _nodes[edge->node(0)];
		for (size_t e = begin; e < end; e++) {
			if (((end - begin) / 3) == 0 || (e - begin) % ((end - begin) / 3) == 0) {
				_corners.push_back(node);
			}
			next(edge, node);
		}
	};

	auto setCorners = [&] (size_t begin, size_t end) {
		if (number == 0 || (end - begin) < 2) {
			return;
		}
		size_t size = end - begin, n = 0, c = number + 1;
		while (n == 0) {
			n = std::ceil((double)size / c--);
		}
		Element *edge = edges[begin], *node = _nodes[edge->node(0)];
		Element *begin_edge = edges[begin], *begin_node = _nodes[edge->node(0)];
		while (edge->domains() == begin_edge->domains() && begin_node->parentEdges().size() == 2) { // find first node
			edge = begin_edge;
			node = begin_node;
			next(begin_edge, begin_node);
		}

		if (begin_node->parentEdges().size() != 2) {
			edge = begin_edge;
			node = begin_node;
			node = node == _nodes[edge->node(0)] ? _nodes[edge->node(edge->nodes() - 1)] : _nodes[edge->node(0)];
		}
		for (size_t i = 1; i < size - 1; i++) {
			if (i % n == 0) {
				_corners.push_back(node);
			}
			next(edge, node);
		}
	};

	size_t begin = 0;
	std::vector<eslocal> edge_nodes = { edges[0]->node(0), edges[0]->node(edges[0]->nodes() - 1) };
	for (size_t i = 1; i < edges.size(); i++) {
		if (edges[i]->domains() != edges[i - 1]->domains()) {
			is_cycle(edge_nodes) ? fixCycle(begin, i) : setCorners(begin, i);
			begin = i;
			edge_nodes.clear();
		}

		if (_nodes[edges[i]->node(0)]->domains().size() > _nodes[edges[i]->node(edges[i]->nodes() - 1)]->domains().size()) {
			_corners.push_back(_nodes[edges[i]->node(0)]);
		}
		if (_nodes[edges[i]->node(0)]->domains().size() < _nodes[edges[i]->node(edges[i]->nodes() - 1)]->domains().size()) {
			_corners.push_back(_nodes[edges[i]->node(edges[i]->nodes() - 1)]);
		}
		edge_nodes.push_back(edges[i]->node(0));
		edge_nodes.push_back(edges[i]->node(edges[i]->nodes() - 1));
	}
	is_cycle(edge_nodes) ? fixCycle(begin, edges.size()) : setCorners(begin, edges.size());

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


void computeDOFsCounters(std::vector<Element*> &elements, const std::vector<Property> &DOFs, const Mesh &mesh)
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


