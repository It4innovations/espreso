#include "mesh.h"
#include "esoutput.h"

using namespace espreso;


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
	std::vector<eslocal> fixPoints(number);
	_fixPoints.resize(parts());

	cilk_for (size_t part = 0; part < parts(); part++) {
		size_t max = (_partPtrs[part + 1] - _partPtrs[part]) / 20 + 1;
		eslocal *eSubPartition = getPartition(_partPtrs[part], _partPtrs[part + 1], std::min(number, max));

		for (eslocal j = 0; j < number; j++) {
			fixPoints[j] = getCentralNode(_partPtrs[part], _partPtrs[part + 1], eSubPartition, part, j);
		}
		std::sort(fixPoints.begin(), fixPoints.end());

		// Remove the same points
		Esutils::unique(fixPoints);

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

void Mesh::saveNodeArray(eslocal *nodeArray, size_t part) const
{
	for (eslocal i = _partPtrs[part]; i < _partPtrs[part + 1]; i++) {
		const Element* e = _elements[i];
		memcpy(nodeArray + i * e->nodes(), e->indices(), e->nodes() * sizeof(eslocal));
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
	ESINFO(GLOBAL_ERROR) << "Check boundaries for surface";
//	surface.computeBoundaries();
//	surface._clusterBoundaries = _clusterBoundaries;
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


static std::vector<eslocal> getIntersection(
		const Element *face,
		const std::vector<std::vector<eslocal> > &boundaries)
{
	std::vector<eslocal> intersection(boundaries[face->node(face->nodes() - 1)]);
	std::vector<eslocal>::iterator it = intersection.end();

	// compute intersection of all nodes
	for (size_t n = face->nodes() - 2; it - intersection.begin() >= 1 &&  n < face->nodes(); n--) {
		std::vector<eslocal> tmp(intersection.begin(), it);
		it = std::set_intersection(tmp.begin(), tmp.end(),
				boundaries[face->node(n)].begin(), boundaries[face->node(n)].end(),
				intersection.begin());
	}

	intersection.resize(it - intersection.begin());
	return intersection;
}

std::vector<std::vector<eslocal> > Mesh::subdomainsInterfaces(Mesh &interface) const
{
	std::vector<std::vector<eslocal> > sMap;

	// 1st step -> create list of faces for each sub-domain

	std::vector<std::vector<Element*> > partFaces(parts());
	std::vector<std::vector<std::vector<eslocal> > > subdomains(parts());
	std::vector<std::vector<std::vector<eslocal> > > elements(parts());

	cilk_for (size_t p = 0; p < parts(); p++) {
		elements[p].resize(_coordinates.localSize(p));
		for (eslocal e = _partPtrs[p]; e < _partPtrs[p + 1]; e++) {
			for (size_t n = 0; n < _elements[e]->nodes(); n++) {
				elements[p][_elements[e]->node(n)].push_back(e);
			}
		}
	}

	cilk_for (size_t p = 0; p < parts(); p++) {
		for (eslocal e = _partPtrs[p]; e < _partPtrs[p + 1]; e++) {
			for (size_t f = 0; f < _elements[e]->faces(); f++) {
				Element* face = _elements[e]->face(f);
				ESINFO(GLOBAL_ERROR) << "subdomainInterfaces";
				std::vector<eslocal> intersection; // = getIntersection(face, _subdomainBoundaries.boundary());
				if (intersection.size() == 1 || intersection[0] != p) {
					delete face; // inner face
				} else {
					bool pass = false;
					std::vector<eslocal> clusterIndices(face->indices(), face->indices() + face->nodes());
					for (size_t i = 1; i < intersection.size(); i++) {
						for (size_t j = 0; j < face->nodes(); j++) {
							face->node(j) = _coordinates.localIndex(clusterIndices[j], intersection[i]);
						}
						if (getIntersection(face, elements[intersection[i]]).size()) {
							if (!pass) {
								partFaces[p].push_back(face);
								subdomains[p].push_back(std::vector<eslocal>(1, intersection[0]));
								pass = true;
							}
							subdomains[p].back().push_back(intersection[i]);
						}
					}
					if (pass) {
						for (size_t j = 0; j < face->nodes(); j++) {
							face->node(j) = clusterIndices[j];
						}
					} else {
						delete face;
					}
				}
			}
		}
	}

	// 2nd step -> make the final mesh
	size_t threads = config::env::CILK_NWORKERS;
	std::vector<size_t> distribution = Esutils::getDistribution(threads, _coordinates.clusterSize());
	std::vector<size_t> offsets(threads);
	std::vector<std::vector<eslocal> > indices(threads);

	cilk_for (size_t t = 0; t < threads; t++) {
		size_t offset = 0;
		indices[t].reserve(distribution[t + 1] - distribution[t]);
		for (size_t i = distribution[t]; i < distribution[t + 1]; i++) {
			indices[t].push_back(_nodes[i]->domains().size() > 1 ? offset++ : -1);
		}
		offsets[t] = offset;
	}

	size_t sum = 0, prev;
	for (size_t t = 0; t < threads; t++) {
		prev = offsets[t];
		offsets[t] = sum;
		sum += prev;
	}

	cilk_for (size_t p = 0; p < parts(); p++) {
		for (size_t i = 0; i < partFaces[p].size(); i++) {
			for (size_t n = 0; n < partFaces[p][i]->nodes(); n++) {
				eslocal offset = partFaces[p][i]->node(n) / distribution[1];
				partFaces[p][i]->node(n) = indices[offset][partFaces[p][i]->node(n) % distribution[1]] + offsets[offset];
			}
		}
	}

	interface._coordinates.reserve(sum);
	size_t index = 0;
	for (size_t t = 0; t < threads; t++) {
		for (size_t i = 0; i < indices[t].size(); i++) {
			if (indices[t][i] >= 0) {
				interface._coordinates.add(_coordinates[distribution[t] + i], index++, distribution[t] + i);
			}
		}
	}

	size_t fSize = 0;
	for (size_t p = 0; p < parts(); p++) {
		fSize += partFaces[p].size();
	}
	interface._elements.reserve(fSize);

	std::vector<std::vector<eslocal> > permutation(parts());
	cilk_for (size_t p = 0; p < parts(); p++) {
		permutation[p].reserve(partFaces[p].size());
		for (size_t i = 0; i < partFaces[p].size(); i++) {
			permutation[p].push_back(i);
		}
		std::sort(permutation[p].begin(), permutation[p].end(), [&] (const eslocal &p1, const eslocal &p2) {
			return subdomains[p][p1] < subdomains[p][p2];
		});
	}

	interface._partPtrs.clear();
	interface._partPtrs.push_back(0);
	for (size_t p = 0; p < parts(); p++) {
		for (size_t i = 0; i < permutation[p].size(); i++) {
			interface._elements.push_back(partFaces[p][permutation[p][i]]);
			if (i + 1 == permutation[p].size() || subdomains[p][permutation[p][i]][1] != subdomains[p][permutation[p][i + 1]][1]) {
				interface._partPtrs.push_back(interface._elements.size());
				interface.makePartContinuous(interface.parts() - 1);
				sMap.insert(sMap.end(), interface.parts() - sMap.size(), subdomains[p][permutation[p][i]]);
			}
		}
	}

	return sMap;
}

void Mesh::computeBorderLinesAndVertices(const Mesh &faces,std::vector<bool> &border, Mesh &lines, std::set<eslocal> &vertices)
{
	lines._elements.clear();
	lines._partPtrs.clear();
	lines._coordinates.clear();
	border.clear();
	border.resize(faces.coordinates().clusterSize(), false);

	std::vector<std::vector<std::vector<eslocal> > > nodesFaces(faces.parts());
	for (size_t p = 0; p < faces.parts(); p++) {
		nodesFaces[p].resize(faces.coordinates().localSize(p));
		for (eslocal e = faces._partPtrs[p]; e < faces._partPtrs[p + 1]; e++) {
			for (size_t n = 0; n < faces._elements[e]->nodes(); n++) {
				// add node's adjacent element
				nodesFaces[p][faces._elements[e]->node(n)].push_back(e);
			}
		}
	}

	std::vector<std::tuple<eslocal, eslocal, eslocal> > commonLines;
	std::vector<char> nSubdomains(faces.coordinates().clusterSize() * faces.parts(), 0);
	std::vector<eslocal> points(faces.coordinates().clusterSize(), 0);
	lines._partPtrs.clear();
	lines._partPtrs.push_back(0);
	for (size_t p = 0; p < faces.parts(); p++) {
		for (size_t e = faces._partPtrs[p]; e < faces._partPtrs[p + 1]; e++) {
			for (size_t f = 0; f < faces._elements[e]->faces(); f++) {
				std::vector<eslocal> line(faces._elements[e]->face(f)->indices(), faces._elements[e]->face(f)->indices() + faces._elements[e]->face(f)->nodes());
				if (isOuterFace(nodesFaces[p], line)) {
					for (size_t n = 0; n < line.size(); n++) {
						line[n] = faces.coordinates().clusterIndex(line[n], p);
						nSubdomains[line[n] * faces.parts() + p] = 1;
						points[line[n]] = 1;
						border[line[n]] = true;
					}
					if (line.front() > line.back()) {
						eslocal tmp = line.front();
						line.front() = line.back();
						line.back() = tmp;
					}
					switch (line.size()) {
					case 2:
						commonLines.push_back(std::make_tuple(line[0], line[1], -1));
						break;
					case 3:
						commonLines.push_back(std::make_tuple(line[0], line[1], line[2]));
						break;
					default:
						ESINFO(ERROR) << "MESH ERROR: unknown line type.";
						exit(EXIT_FAILURE);
					}
				}
			}
		}
	}
	std::sort(commonLines.begin(), commonLines.end());
	auto it = std::unique(commonLines.begin(), commonLines.end());
	commonLines.resize(std::distance(commonLines.begin(), it));

	eslocal linePoints = 0;
	for (size_t i = 0; i < faces.coordinates().clusterSize(); i++) {
		if (points[i] == 1) {
			lines._coordinates.add(faces.coordinates()[i], linePoints, i);
			points[i] = linePoints++;
		}
	}

	// commonLines -> unique common lines in global indexing
	// points -> mapping from global to cluster indices
	// nSubdomains -> map: global node -> parts


	std::vector<bool> procesedPoints(faces.coordinates().clusterSize(), false);
	size_t processedCounter = 0;
	// skip all points not on border
	for (size_t i = 0; i < procesedPoints.size(); i++) {
		char *s = nSubdomains.data() + i * faces.parts();
		char *e = nSubdomains.data() + (i + 1) * faces.parts();
		if (std::all_of(s, e, [](char flag) { return flag == 0; })) {
			procesedPoints[i] = true;
			processedCounter++;
		}
	}

	// check whether two points are in the same sub-domains
	auto same_subdomains = [&] (size_t i, size_t j) {
		char *s = nSubdomains.data() + i * faces.parts();
		char *e = nSubdomains.data() + (i + 1) * faces.parts();
		char *t = nSubdomains.data() + j * faces.parts();
		return std::equal(s, e, t);
	};

	auto subdomains_count = [&] (eslocal i) {
		return std::count(nSubdomains.data() + i * faces.parts(), nSubdomains.data() + (i + 1) * faces.parts(), 1);
	};

	size_t begin = 0;
	std::vector<bool> cyclic_line;
	while (processedCounter < procesedPoints.size()) {
		// skip to the first unprocessed point
		while (procesedPoints[begin]) {
			begin++;
		}

		// process all points on actual border
		for (size_t i = begin; i < procesedPoints.size(); i++) {
			if (same_subdomains(begin, i)) {
				procesedPoints[i] = true;
				processedCounter++;
			}
		}

		bool innerLine = subdomains_count(begin) > 1;
		bool flag = true;
		eslocal params[6] = {0, 0, 0, 0, 0, 0};
		for (size_t i = 0; i < commonLines.size(); i++) {
			eslocal start = std::get<0>(commonLines[i]);
			eslocal mid = std::get<1>(commonLines[i]);
			eslocal end = (std::get<2>(commonLines[i]) == -1) ? mid : std::get<2>(commonLines[i]);

			if (same_subdomains(begin, start) && same_subdomains(begin, end)) {
				eslocal tmp[3] = { points[start], points[mid], points[end] };
				if (std::get<2>(commonLines[i]) == -1) {
					lines._elements.push_back(new Line2(tmp));
				} else {
					lines._elements.push_back(new Line3(tmp));
				}
			}
			if (same_subdomains(begin, start) != same_subdomains(begin, end)) {
				if (subdomains_count(start) > subdomains_count(end)) {
					vertices.insert(points[start]);
				} else {
					vertices.insert(points[end]);
				}
			}
		}
		for (eslocal e = lines._elements.size() - 1; e >= lines._partPtrs.back(); e--) {
			if (vertices.find(lines._elements[e]->node(0)) != vertices.end() ||
				vertices.find(lines._elements[e]->node(lines._elements[e]->nodes() - 1)) != vertices.end()) {

				// erase elements with vertex corners
				lines._elements.erase(lines._elements.begin() + e);
			}
		}

		if (lines._partPtrs.back() < lines._elements.size()) {
			lines._partPtrs.push_back(lines._elements.size());
			lines.makePartContinuous(lines.parts() - 1);
		}
	}

	lines.mapCoordinatesToDomains();
}

void Mesh::correctCycle(Mesh &faces, Mesh &lines, bool average)
{
	auto checkCycle = [&] (eslocal part) {
		std::vector<int> counter(lines._coordinates.clusterSize());
		for (size_t e = lines._partPtrs[part]; e < lines._partPtrs[part + 1]; e++) {
			counter[lines._elements[e]->node(0)]++;
			counter[lines._elements[e]->node(lines._elements[e]->nodes() - 1)]++;
		}
		return std::all_of(counter.begin(), counter.end(), [] (int count) { return count % 2 == 0; });
	};

	std::vector<eslocal> partPtrs(lines._partPtrs);
	eslocal offset = 0;
	for (size_t p = 0; p < lines.parts(); p++) {
		if (!checkCycle(p)) {
			continue;
		}

		eslocal begin = lines._partPtrs[p], end = lines._partPtrs[p + 1];

		// defense strategy -> it is better to add more corners than is necessary
		int max = end - begin;
		size_t corners = std::min(max, 4);
		eslocal *ePartition = lines.getPartition(begin, end, corners);

		if (average) {
			// divide to sub-parts -> corners will be add later
			eslocal counter = 0;
			for (size_t i = 0; i < corners; i++) {
				for (size_t e = 0; e < end - begin; e++) {
					if (ePartition[e] == i) {
						std::swap(ePartition[e], ePartition[counter]);
						std::swap(lines._elements[begin + e], lines._elements[begin + counter]);
						counter++;
					}
				}
				if (i + 1 < corners) {
					partPtrs.insert(partPtrs.begin() + p + 1 + offset, begin + counter);
					offset++;
				}
			}
		} else {
			// set corners directly and remove line
			for (eslocal j = 0; j < corners; j++) {
				eslocal sCorner = lines.getCentralNode(begin, end, ePartition, p, j);
				sCorner = lines.coordinates().globalIndex(sCorner, p);
				sCorner = faces.coordinates().globalIndex(sCorner);
				ESINFO(GLOBAL_ERROR) << "set corner";
				//_subdomainBoundaries.setCorner(sCorner);
			}
			for (size_t e = begin; e < end; e++) {
				delete lines._elements[e];
				lines._elements[e] = NULL;
			}

			for (size_t pp = p + 1 + offset; pp < lines.parts() + offset; pp++) {
				partPtrs[pp] = partPtrs[pp + 1] - (end - begin);
			}
			partPtrs.pop_back();
			offset--;
		}
		delete[] ePartition;
	}

	if (!average) {
		for (eslocal p = lines.parts() - 1; p >= 0; p--) {
			if (lines._elements[lines._partPtrs[p]] == NULL) {
				lines._coordinates._clusterIndex.erase(lines._coordinates._clusterIndex.begin() + p);
				lines._elements.erase(lines._elements.begin() + lines._partPtrs[p], lines._elements.begin() + lines._partPtrs[p + 1]);
			}
		}
	}
	lines._partPtrs.swap(partPtrs);
	if (average) {
		lines.mapCoordinatesToDomains();
		//lines.computeFixPoints(1);
	}
}

void Mesh::prepareAveragingLines(Mesh &faces, Mesh &lines)
{
	// check whether a point has dirichler condition
	auto has_dirichlet = [&] (eslocal i) {
		return false;
//		eslocal index = lines.coordinates().globalIndex(i);
//		index = faces.coordinates().globalIndex(index);
//		auto &dx = coordinates().property(DIRICHLET_X).values();
//		auto &dy = coordinates().property(DIRICHLET_Y).values();
//		auto &dz = coordinates().property(DIRICHLET_Z).values();
//		return dx.find(index) != dx.end() || dy.find(index) != dy.end() || dz.find(index) != dz.end();
	};

	// check whether a point is on cluster boundary
	auto on_cluster_boundary = [&] (eslocal i) {
		eslocal index = lines.coordinates().globalIndex(i);
		index = faces.coordinates().globalIndex(index);
		return _nodes[index]->clusters().size() > 1;
	};

	//lines.remapElementsToCluster();

	for (size_t e = 0; e < lines._partPtrs.back(); e++) {
		for (size_t n = 0; n < lines._elements[e]->nodes(); n++) {
			if (has_dirichlet(lines._elements[e]->node(n)) || on_cluster_boundary(lines._elements[e]->node(n))) {
				delete lines._elements[e];
				lines._elements[e] = NULL;
				break;
			}
		}
	}

	size_t counter = 0;
	std::vector<eslocal> partition(1, 0);
	for (size_t p = 0; p < lines.parts(); p++) {
		for (size_t e = lines._partPtrs[p]; e < lines._partPtrs[p + 1]; e++) {
			if (lines._elements[e] != NULL) {
				std::swap(lines._elements[counter++], lines._elements[e]);
			}
		}
		if (counter > partition.back()) {
			partition.push_back(counter);
		}
	}
	lines._elements.resize(partition.back());
	lines._partPtrs.swap(partition);

	lines.mapCoordinatesToDomains();
}

void Mesh::prepareAveragingFaces(Mesh &faces, std::vector<bool> &border)
{
	// faces.remapElementsToCluster();

	for (size_t e = 0; e < faces._partPtrs.back(); e++) {
		for (size_t n = 0; n < faces._elements[e]->nodes(); n++) {
			if (border[faces._elements[e]->node(n)]) {
				delete faces._elements[e];
				faces._elements[e] = NULL;
				break;
			}
		}
	}

	size_t counter = 0;
	std::vector<eslocal> partition(1, 0);
	for (size_t p = 0; p < faces.parts(); p++) {
		for (size_t e = faces._partPtrs[p]; e < faces._partPtrs[p + 1]; e++) {
			if (faces._elements[e] != NULL) {
				std::swap(faces._elements[counter++], faces._elements[e]);
			}
		}
		if (counter > partition.back()) {
			partition.push_back(counter);
		}
	}
	faces._elements.resize(partition.back());
	faces._partPtrs.swap(partition);

	faces.mapCoordinatesToDomains();
}


void Mesh::computeCorners(eslocal number, bool vertices, bool edges, bool faces)
{
	if (parts() == 1 || (!vertices && !edges && !faces)) {
		return;
	}

	computeFacesSharedByDomains();
	computeEdgesOnBordersOfFacesSharedByDomains();


//	if (parts() < 1) {
//		ESINFO(ERROR) << "Internal error: _partPtrs.size().";
//		exit(EXIT_FAILURE);
//	}
//	if (parts() == 1 || (!vertices && !edges && !faces && !averageEdges && !averageFaces)) {
//		return;
//	}
//
//	_subdomainBoundaries._corners.clear();
//	_subdomainBoundaries._corners.resize(_subdomainBoundaries.size(), false);
//
//	Mesh commonFaces;
//	Mesh commonLines;
//	std::set<eslocal> commonVertices;
//	std::vector<bool> commonFacesBorder;
//
//	subdomainsInterfaces(commonFaces);
//	computeBorderLinesAndVertices(commonFaces, commonFacesBorder, commonLines, commonVertices);
//
//	if (config::output::SAVE_FACES) {
//		output::VTK_Full::mesh(commonFaces, "meshFaces", config::output::SUBDOMAINS_SHRINK_RATIO, config::output::CLUSTERS_SHRINK_RATIO);
//	}
//
//	if (config::output::SAVE_LINES) {
//		output::VTK_Full::mesh(commonLines, "meshLines", config::output::SUBDOMAINS_SHRINK_RATIO, config::output::CLUSTERS_SHRINK_RATIO);
//	}
//
//	auto faceToCluster = [&] (eslocal index, eslocal part) {
//		return commonFaces.coordinates().globalIndex(index, part);
//	};
//
//	auto lineToCluster = [&] (eslocal index, eslocal part) {
//		index = commonLines.coordinates().globalIndex(index, part);
//		return commonFaces.coordinates().globalIndex(index);
//	};
//
//	auto corners = [&] (Mesh &mesh, std::function<eslocal(eslocal, eslocal)> map) {
//		for (size_t p = 0; p < mesh.parts(); p++) {
//			std::vector<eslocal> fixPoints = mesh.computeFixPoints(p, number);
//			for (size_t i = 0; i < fixPoints.size(); i++) {
//				_subdomainBoundaries.setCorner(map(fixPoints[i], p));
//			}
//		}
//	};
//
//	auto average = [&] (Mesh &mesh, std::function<eslocal(eslocal, eslocal)> map) {
//		for (size_t p = 0; p < mesh.parts(); p++) {
//			std::vector<eslocal> fixPoints = mesh.computeFixPoints(p, 1);
//			eslocal corner = map(fixPoints[0], p);
//			_subdomainBoundaries.setCorner(corner);
//
//			std::set<eslocal> aPoints;
//			for (size_t e = mesh._partPtrs[p]; e < mesh._partPtrs[p + 1]; e++) {
//				for (size_t n = 0; n < mesh._elements[e]->size(); n++) {
//					eslocal point = map(mesh._elements[e]->node(n), p);
//					if (!_subdomainBoundaries.isCorner(point)) {
//						aPoints.insert(point);
//					}
//				}
//			}
//
//			std::vector<eslocal>& averaging = _subdomainBoundaries.averaging(corner);
//			averaging.insert(averaging.begin(), aPoints.begin(), aPoints.end());
//		}
//	};
//
//	if (vertices) {
//		for (auto it = commonVertices.begin(); it != commonVertices.end(); ++it) {
//			_subdomainBoundaries.setCorner(commonFaces.coordinates().globalIndex(commonLines.coordinates().globalIndex(*it)));
//		}
//		correctCycle(commonFaces, commonLines, averageEdges);
//		ESINFO(DETAILS) << "Set corners to vertices";
//	}
//
//	if (edges && !averageEdges) {
//		corners(commonLines, lineToCluster);
//		ESINFO(DETAILS) << "Set corners to edges - " << number << " per edge";
//	}
//
//	if (averageEdges) {
//		prepareAveragingLines(commonFaces, commonLines);
//		average(commonLines, lineToCluster);
//		ESINFO(DETAILS) << "Average edged";
//	}
//
//	if (faces) {
//		corners(commonFaces, faceToCluster);
//		ESINFO(DETAILS) << "Set corners to faces - " << number << " per face";
//	}
//	if (averageFaces) {
//		prepareAveragingFaces(commonFaces, commonFacesBorder);
//		average(commonFaces, faceToCluster);
//		ESINFO(DETAILS) << "Average faces";
//	}
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
				if (filter(_nodes, edge)) {
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
				if (filter(_nodes, face)) {
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

	if (config::output::SAVE_FACES) {
		Mesh mesh;
		mesh._coordinates = _coordinates;
		for (size_t i = 0; i < _faces.size(); i++) {
			mesh._elements.push_back(_faces[i]->copy());
		}
		mesh._partPtrs.clear();
		mesh._partPtrs.push_back(0);
		std::sort(mesh._elements.begin(), mesh._elements.end(), [] (Element* e1, Element* e2) { return e1->domains() < e2->domains(); });
		for (size_t i = 1; i < mesh._elements.size(); i++) {
			if (mesh._elements[i]->domains() != mesh._elements[i - 1]->domains()) {
				mesh._partPtrs.push_back(i);
			}
		}
		mesh._partPtrs.push_back(mesh._elements.size());
		mesh.mapElementsToDomains();
		mesh.mapCoordinatesToDomains();

		output::VTK_Full::mesh(mesh, "meshFaces", config::output::SUBDOMAINS_SHRINK_RATIO, config::output::CLUSTERS_SHRINK_RATIO);
	}
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
	for (size_t e = 0; e < _elements.size(); e++) {
		for (size_t n = 0; n < _elements[e]->nodes(); n++) {
			_nodes[_elements[e]->node(n)]->parentElements().push_back(_elements[e]);
		}
	}
}

void Mesh::fillParentEdgesToNodes()
{
	for (size_t e = 0; e < _edges.size(); e++) {
		for (size_t n = 0; n < _edges[e]->nodes(); n++) {
			_nodes[_edges[e]->node(n)]->parentElements().push_back(_edges[e]);
		}
	}
}

void Mesh::fillParentFacesToNodes()
{
	for (size_t f = 0; f < _faces.size(); f++) {
		for (size_t n = 0; n < _faces[f]->nodes(); n++) {
			_nodes[_faces[f]->node(n)]->parentFaces().push_back(_faces[f]);
		}
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
				if (filter(_nodes, edge)) {
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
				if (e2->parentElements()[0]->edge(i) != NULL && *(e1) == *(e2->parentElements()[e]->edge(i))) {
					e2->parentElements()[e]->setEdge(i, e1);
					break;
				}
			}
		}
	});

	mapEdgesToClusters();
	mapEdgesToDomains();

	if (config::output::SAVE_EDGES) {
		Mesh mesh;
		mesh._coordinates = _coordinates;
		for (size_t i = 0; i < _edges.size(); i++) {
			mesh._elements.push_back(_edges[i]->copy());
		}
		mesh._partPtrs.clear();
		mesh._partPtrs.push_back(0);
		std::sort(mesh._elements.begin(), mesh._elements.end(), [] (Element* e1, Element* e2) { return e1->domains() < e2->domains(); });
		for (size_t i = 1; i < mesh._elements.size(); i++) {
			if (mesh._elements[i]->domains() != mesh._elements[i - 1]->domains()) {
				mesh._partPtrs.push_back(i);
			}
		}
		mesh._partPtrs.push_back(mesh._elements.size());
		mesh.mapElementsToDomains();
		mesh.mapCoordinatesToDomains();

		output::VTK_Full::mesh(mesh, "meshEdges", config::output::SUBDOMAINS_SHRINK_RATIO, config::output::CLUSTERS_SHRINK_RATIO);
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
		return intersection.size() == 1;
	});
}

void Mesh::computeFacesSharedByDomains()
{
	fillFacesFromElements([] (const std::vector<Element*> &nodes, const Element *face) {
		std::vector<eslocal> intersection(nodes[face->node(face->nodes() - 1)]->domains()); // it is better to start from end (from mid points)
		auto it = intersection.end();

		for (size_t n = face->nodes() - 2; it - intersection.begin() > 1 &&  n < face->nodes(); n--) {
			std::vector<eslocal> tmp(intersection.begin(), it);
			it = std::set_intersection(tmp.begin(), tmp.end(),
					nodes[face->node(n)]->domains().begin(), nodes[face->node(n)]->domains().end(),
					intersection.begin());
		}

		intersection.resize(it - intersection.begin());
		return intersection.size() > 1;
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

void Mesh::computeEdgesOnBordersOfFacesSharedByDomains()
{
	fillParentFacesToNodes();

	fillEdgesFromFaces([] (const std::vector<Element*> &nodes, const Element* edge) {
//		if (edge->parentElements()[0]->domains().size() != 2) {
//			return false;
//		}
//
//		std::vector<eslocal> intersection(nodes[edge->node(edge->nodes() - 1)]->domains()); // it is better to start from end (from mid points)
//		auto it = intersection.end();
//
//		for (size_t n = edge->nodes() - 2; it - intersection.begin() > 2 &&  n < edge->nodes(); n--) {
//			std::vector<eslocal> tmp(intersection.begin(), it);
//			it = std::set_intersection(tmp.begin(), tmp.end(),
//					nodes[edge->node(n)]->domains().begin(), nodes[edge->node(n)]->domains().end(),
//					intersection.begin());
//		}
//
//		intersection.resize(it - intersection.begin());
//		if (intersection.size() > 2) {
//			return true;
//		}

		for (size_t n = edge->nodes() - 1; n < edge->nodes(); n--) {
			if (nodes[edge->node(n)]->parentFaces().size() > 2) {
//				for (size_t f = 1; f < nodes[edge->node(n)]->parentFaces().size(); f++) {
//					if (nodes[edge->node(n)]->parentFaces()[f]->domains() != nodes[edge->node(n)]->parentFaces()[f - 1]->domains()) {
//						return true;
//					}
//				}
				return false;
			}
		}
		return true;
	});

	fillParentEdgesToNodes();
}

void Mesh::clearEdgesWithoutSettings()
{

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
				elements[i]->domains().push_back(elements[i]->parentElements()[e]->domains()[0]);
			}
			std::sort(elements[i]->domains().begin(), elements[i]->domains().end());
			Esutils::unique(elements[i]->domains());
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
	int rank = 0;

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

					sBuffer[t][*c].push_back(elements[e]->vtkCode());
					for (size_t n = 0; n < elements[e]->coarseNodes(); n++) {
						sBuffer[t][*c].push_back(mesh.coordinates().globalIndex(elements[e]->node(n)));
					}

					for (size_t i = 0; i < DOFs.size(); i++) {
						sBuffer[t][*c].push_back(elements[e]->DOFsDomainsCounters()[cluster * DOFs.size() + i]);
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
				size_t cluster = std::lower_bound((*it)->clusters().begin(), (*it)->clusters().end(), n) - (*it)->clusters().begin();
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

	cilk_for (size_t p = 0; p < parts(); p++) {
		std::vector<eslocal> l2g;
		for (eslocal e = _partPtrs[p]; e < _partPtrs[p + 1]; e++) {
			l2g.insert(l2g.end(), _elements[e]->indices(), _elements[e]->indices() + _elements[e]->nodes());
		}

		std::sort(l2g.begin(), l2g.end());
		Esutils::unique(l2g);

		_coordinates._clusterIndex[p] = l2g;
	}
}

}


