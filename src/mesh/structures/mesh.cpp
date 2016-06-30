
#include "mesh.h"
#include "esoutput.h"

namespace espreso {

Mesh::Mesh():_elements(0), _fixPoints(0), _DOFs(3)
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
		remapElementsToSubdomain();
		computeFixPoints(0);
		computeBoundaries();
		return;
	}

	if (this->parts()) {
		this->remapElementsToCluster();
	}

	_partPtrs.resize(parts + 1);
	_partPtrs[0] = 0;

	eslocal *ePartition = getPartition(0, _elements.size(), parts);

	Element *e;
	eslocal p;
	for (size_t part = 0; part < parts; part++) {
		eslocal index = _partPtrs[part];	// index of last ordered element
		for (size_t i = _partPtrs[part]; i < _elements.size(); i++) {
			if (ePartition[i] == part) {
				if (i == index) {
					index++;
				} else {
					e = _elements[i];
					_elements[i] = _elements[index];
					_elements[index] = e;
					p = ePartition[i];
					ePartition[i] = ePartition[index];
					ePartition[index] = p;
					index++;
				}
			}
		}
		_partPtrs[part + 1] = index;
		ESTEST(MANDATORY) << "subdomain without element" << (_partPtrs[part] == _partPtrs[part + 1] ? TEST_FAILED : TEST_PASSED);
	}
	delete[] ePartition;

	remapElementsToSubdomain();
	computeFixPoints(0);
	computeBoundaries();
}

void APIMesh::partitiate(size_t parts)
{
	if (parts == 1 && this->parts() == 1) {
		_partPtrs.resize(parts + 1);
		_partPtrs[0] = 0;
		_partPtrs[1] = _elements.size();
		remapElementsToSubdomain();
		computeFixPoints(0);
		computeBoundaries();
		return;
	}

	if (this->parts()) {
		this->remapElementsToCluster();
	}

	_partPtrs.resize(parts + 1);
	_partPtrs[0] = 0;

	eslocal *ePartition = getPartition(0, _elements.size(), parts);

	Element *e;
	eslocal p;
	for (size_t part = 0; part < parts; part++) {
		eslocal index = _partPtrs[part];	// index of last ordered element
		for (size_t i = _partPtrs[part]; i < _elements.size(); i++) {
			if (ePartition[i] == part) {
				if (i == index) {
					index++;
				} else {
					e = _elements[i];
					_elements[i] = _elements[index];
					_elements[index] = e;
					p = ePartition[i];
					ePartition[i] = ePartition[index];
					ePartition[index] = p;
					index++;
					_eMatrices[i].swap(_eMatrices[index]);
				}
			}
		}
		_partPtrs[part + 1] = index;
		ESTEST(MANDATORY) << "subdomain without element" << (_partPtrs[part] == _partPtrs[part + 1] ? TEST_FAILED : TEST_PASSED);
	}
	delete[] ePartition;

	remapElementsToSubdomain();
	computeFixPoints(0);
	computeBoundaries();
}

void Mesh::computeBoundaries()
{
	_subdomainBoundaries.clear();
	_subdomainBoundaries.resize(_coordinates.clusterSize());

	// reset corners
	_subdomainBoundaries._averaging.clear();
	std::fill(_subdomainBoundaries._corners.begin(), _subdomainBoundaries._corners.end(), false);

	for (size_t p = 0; p < parts(); p++) {
		for (eslocal e = _partPtrs[p]; e < _partPtrs[p + 1]; e++) {
			for (size_t n = 0; n < _elements[e]->size(); n++) {
				auto index = _coordinates.clusterIndex(_elements[e]->node(n), p);
				if (!_subdomainBoundaries[index].size() || _subdomainBoundaries[index].back() != p) {
					_subdomainBoundaries[index].push_back(p);
				}
			}
		}
	}
}

void Mesh::computeFixPoints(size_t fixPoints)
{
	_fixPoints.clear();
	_fixPoints.resize(parts(), std::vector<eslocal>(fixPoints));

	if (fixPoints == 0) {
		return;
	}

	cilk_for (eslocal i = 0; i < parts(); i++) {
		size_t max = (_partPtrs[i + 1] - _partPtrs[i]) / 20 + 1;
		eslocal *eSubPartition = getPartition(_partPtrs[i], _partPtrs[i + 1], std::min(fixPoints, max));

		for (eslocal j = 0; j < fixPoints; j++) {
			_fixPoints[i][j] = getCentralNode(_partPtrs[i], _partPtrs[i + 1], eSubPartition, i, j);
		}
		std::sort(_fixPoints[i].begin(), _fixPoints[i].end());

		// Remove the same points
		auto it = std::unique(_fixPoints[i].begin(), _fixPoints[i].end());
		_fixPoints[i].resize(it - _fixPoints[i].begin());

		delete[] eSubPartition;
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
		e[index + 1] = e[index] + _elements[i]->coarseSize();
	}

	// create array of nodes
	n = new eslocal[e[eSize]];
	// number of common nodes to be neighbor
	ncommon = 4;
	for (eslocal i = first, index = 0; i < last; i++, index++) {
		const Element* el = _elements[i];
		memcpy(n + e[index], el->indices(), el->coarseSize() * sizeof(eslocal));
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
			for (size_t j = 0; j < _elements[i]->size(); j++) {
				std::vector<eslocal> neigh = _elements[i]->getNeighbours(j);
				for (size_t k = 0; k < neigh.size(); k++) {
					if (_elements[i]->node(j) < neigh[k]) {
						neighbours[_elements[i]->node(j)].insert(neigh[k]);
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

	return result;
}

Mesh::~Mesh()
{
	for (size_t i = 0; i < _elements.size(); i++) {
		delete _elements[i];
	}
}

void Mesh::saveNodeArray(eslocal *nodeArray, size_t part) const
{
	for (eslocal i = _partPtrs[part]; i < _partPtrs[part + 1]; i++) {
		const Element* e = _elements[i];
		memcpy(nodeArray + i * e->size(), e->indices(), e->size() * sizeof(eslocal));
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
			for (size_t k = 0; k < _elements[j]->size(); k++) {
				nodesElements[_elements[j]->node(k)].push_back(j);
			}
		}

		// compute number of elements and fill used nodes
		for (eslocal j = _partPtrs[i]; j < _partPtrs[i + 1]; j++) {
			for (size_t k = 0; k < _elements[j]->faces(); k++) {
				std::vector<eslocal> face = _elements[j]->getFace(k);
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

	surface.coordinates() = _coordinates;

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

	surface.remapElementsToSubdomain();
	surface.computeFixPoints(0);
	surface.computeBoundaries();
	surface._clusterBoundaries = _clusterBoundaries;
	surface._DOFs = _DOFs;
	surface._neighbours = _neighbours;
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
		ePtr[i] = ePtr[i - 1] + _elements[e]->coarseSize();
	}

	nCommon = 4;
	eInd = new eslocal[ePtr[ne]];
	for (size_t e = begin, i = 0; e < end; e++, i++) {
		const Element* el = _elements[e];
		memcpy(eInd + ePtr[i], el->indices(), el->coarseSize() * sizeof(eslocal));
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
	std::vector<eslocal> intersection(boundaries[face->node(face->size() - 1)]);
	std::vector<eslocal>::iterator it = intersection.end();

	// compute intersection of all nodes
	for (size_t n = face->size() - 2; it - intersection.begin() >= 1 &&  n < face->size(); n--) {
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
			for (size_t n = 0; n < _elements[e]->size(); n++) {
				elements[p][_elements[e]->node(n)].push_back(e);
			}
		}
	}

	remapElementsToCluster();

	cilk_for (size_t p = 0; p < parts(); p++) {
		for (eslocal e = _partPtrs[p]; e < _partPtrs[p + 1]; e++) {
			for (size_t f = 0; f < _elements[e]->faces(); f++) {
				Element* face = _elements[e]->getFullFace(f);
				std::vector<eslocal> intersection = getIntersection(face, _subdomainBoundaries.boundary());
				if (intersection.size() == 1 || intersection[0] != p) {
					delete face; // inner face
				} else {
					bool pass = false;
					std::vector<eslocal> clusterIndices(face->indices(), face->indices() + face->size());
					for (size_t i = 1; i < intersection.size(); i++) {
						for (size_t j = 0; j < face->size(); j++) {
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
						for (size_t j = 0; j < face->size(); j++) {
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
			indices[t].push_back(_subdomainBoundaries[i].size() > 1 ? offset++ : -1);
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
			for (size_t n = 0; n < partFaces[p][i]->size(); n++) {
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

	interface.remapElementsToSubdomain();
	interface.computeFixPoints(0);
	interface.computeBoundaries();

	remapElementsToSubdomain();

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
			for (size_t n = 0; n < faces._elements[e]->size(); n++) {
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
				std::vector<eslocal> line = faces._elements[e]->getFace(f);
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
					lines._elements.push_back(new Line2(tmp, params));
				} else {
					lines._elements.push_back(new Line3(tmp, params));
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
				vertices.find(lines._elements[e]->node(lines._elements[e]->size() - 1)) != vertices.end()) {

				// erase elements with vertex corners
				lines._elements.erase(lines._elements.begin() + e);
			}
		}

		if (lines._partPtrs.back() < lines._elements.size()) {
			lines._partPtrs.push_back(lines._elements.size());
			lines.makePartContinuous(lines.parts() - 1);
		}
	}

	lines.remapElementsToSubdomain();
	lines.computeFixPoints(0);
}

void Mesh::correctCycle(Mesh &faces, Mesh &lines, bool average)
{
	auto checkCycle = [&] (eslocal part) {
		std::vector<int> counter(lines._coordinates.clusterSize());
		for (size_t e = lines._partPtrs[part]; e < lines._partPtrs[part + 1]; e++) {
			counter[lines._elements[e]->node(0)]++;
			counter[lines._elements[e]->node(lines._elements[e]->size() - 1)]++;
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
				_subdomainBoundaries.setCorner(sCorner);
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
	} else {
		lines.remapElementsToCluster();
	}
	lines._partPtrs.swap(partPtrs);
	if (average) {
		lines.remapElementsToSubdomain();
		lines.computeFixPoints(1);
	}
}

void Mesh::prepareAveragingLines(Mesh &faces, Mesh &lines)
{
	// check whether a point has dirichler condition
	auto has_dirichlet = [&] (eslocal i) {
		eslocal index = lines.coordinates().globalIndex(i);
		index = faces.coordinates().globalIndex(index);
		auto &dx = coordinates().property(DIRICHLET_X).values();
		auto &dy = coordinates().property(DIRICHLET_Y).values();
		auto &dz = coordinates().property(DIRICHLET_Z).values();
		return dx.find(index) != dx.end() || dy.find(index) != dy.end() || dz.find(index) != dz.end();
	};

	// check whether a point is on cluster boundary
	auto on_cluster_boundary = [&] (eslocal i) {
		eslocal index = lines.coordinates().globalIndex(i);
		index = faces.coordinates().globalIndex(index);
		auto &cb = clusterBoundaries();
		return cb[index].size() > 1;
	};

	lines.remapElementsToCluster();

	for (size_t e = 0; e < lines._partPtrs.back(); e++) {
		for (size_t n = 0; n < lines._elements[e]->size(); n++) {
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

	lines.remapElementsToSubdomain();
}

void Mesh::prepareAveragingFaces(Mesh &faces, std::vector<bool> &border)
{
	faces.remapElementsToCluster();

	for (size_t e = 0; e < faces._partPtrs.back(); e++) {
		for (size_t n = 0; n < faces._elements[e]->size(); n++) {
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

	faces.remapElementsToSubdomain();
}


void Mesh::computeCorners(eslocal number, bool vertices, bool edges, bool faces, bool averageEdges, bool averageFaces)
{
	if (parts() < 1) {
		ESINFO(ERROR) << "Internal error: _partPtrs.size().";
		exit(EXIT_FAILURE);
	}
	if (parts() == 1 || (!vertices && !edges && !faces && !averageEdges && !averageFaces)) {
		return;
	}

	_subdomainBoundaries._corners.clear();
	_subdomainBoundaries._corners.resize(_subdomainBoundaries.size(), false);

	Mesh commonFaces;
	Mesh commonLines;
	std::set<eslocal> commonVertices;
	std::vector<bool> commonFacesBorder;

	subdomainsInterfaces(commonFaces);
	computeBorderLinesAndVertices(commonFaces, commonFacesBorder, commonLines, commonVertices);

	if (config::output::SAVE_FACES) {
		output::VTK_Full::mesh(commonFaces, "meshFaces", config::output::SUBDOMAINS_SHRINK_RATIO, config::output::CLUSTERS_SHRINK_RATIO);
	}

	if (config::output::SAVE_LINES) {
		output::VTK_Full::mesh(commonLines, "meshLines", config::output::SUBDOMAINS_SHRINK_RATIO, config::output::CLUSTERS_SHRINK_RATIO);
	}

	auto faceToCluster = [&] (eslocal index, eslocal part) {
		return commonFaces.coordinates().globalIndex(index, part);
	};

	auto lineToCluster = [&] (eslocal index, eslocal part) {
		index = commonLines.coordinates().globalIndex(index, part);
		return commonFaces.coordinates().globalIndex(index);
	};

	auto corners = [&] (Mesh &mesh, std::function<eslocal(eslocal, eslocal)> map) {
		mesh.computeFixPoints(number);
		for (size_t p = 0; p < mesh.parts(); p++) {
			for (size_t i = 0; i < mesh.getFixPoints()[p].size(); i++) {
				_subdomainBoundaries.setCorner(map(mesh.getFixPoints()[p][i], p));
			}
		}
	};

	auto average = [&] (Mesh &mesh, std::function<eslocal(eslocal, eslocal)> map) {
		mesh.computeFixPoints(1);
		for (size_t p = 0; p < mesh.parts(); p++) {
			eslocal corner = map(mesh._fixPoints[p][0], p);
			_subdomainBoundaries.setCorner(corner);

			std::set<eslocal> aPoints;
			for (size_t e = mesh._partPtrs[p]; e < mesh._partPtrs[p + 1]; e++) {
				for (size_t n = 0; n < mesh._elements[e]->size(); n++) {
					eslocal point = map(mesh._elements[e]->node(n), p);
					if (!_subdomainBoundaries.isCorner(point)) {
						aPoints.insert(point);
					}
				}
			}

			std::vector<eslocal>& averaging = _subdomainBoundaries.averaging(corner);
			averaging.insert(averaging.begin(), aPoints.begin(), aPoints.end());
		}
	};

	if (vertices) {
		for (auto it = commonVertices.begin(); it != commonVertices.end(); ++it) {
			_subdomainBoundaries.setCorner(commonFaces.coordinates().globalIndex(commonLines.coordinates().globalIndex(*it)));
		}
		correctCycle(commonFaces, commonLines, averageEdges);
		ESINFO(DETAILS) << "Set corners to vertices";
	}

	if (edges && !averageEdges) {
		corners(commonLines, lineToCluster);
		ESINFO(DETAILS) << "Set corners to edges - " << number << " per edge";
	}

	if (averageEdges) {
		prepareAveragingLines(commonFaces, commonLines);
		average(commonLines, lineToCluster);
		ESINFO(DETAILS) << "Average edged";
	}

	if (faces) {
		corners(commonFaces, faceToCluster);
		ESINFO(DETAILS) << "Set corners to faces - " << number << " per face";
	}
	if (averageFaces) {
		prepareAveragingFaces(commonFaces, commonFacesBorder);
		average(commonFaces, faceToCluster);
		ESINFO(DETAILS) << "Average faces";
	}
}

void Mesh::remapElementsToSubdomain() const
{
	_coordinates._clusterIndex.clear();
	_coordinates._clusterIndex.resize(parts());

	cilk_for (size_t p = 0; p < parts(); p++) {
		std::vector<eslocal> l2g;

		for (eslocal e = _partPtrs[p]; e < _partPtrs[p + 1]; e++) {
			l2g.insert(l2g.end(), _elements[e]->indices(), _elements[e]->indices() + _elements[e]->size());
		}

		std::sort(l2g.begin(), l2g.end());
		auto it = std::unique(l2g.begin(), l2g.end());
		l2g.resize(std::distance(l2g.begin(), it));

		for (eslocal e = _partPtrs[p]; e < _partPtrs[p + 1]; e++) {
			for (eslocal n = 0; n < _elements[e]->size(); n++) {
				_elements[e]->node(n) = std::lower_bound(l2g.begin(), l2g.end(), _elements[e]->node(n)) - l2g.begin();
			}
		}

		_coordinates._clusterIndex[p] = l2g;
	}
}

void Mesh::remapElementsToCluster() const
{
	cilk_for (eslocal p = 0; p < parts(); p++) {
		for (eslocal e = _partPtrs[p]; e < _partPtrs[p + 1]; e++) {
			for (eslocal n = 0; n < _elements[e]->size(); n++) {
				_elements[e]->node(n) = _coordinates.clusterIndex(_elements[e]->node(n), p);
			}
		}
	}
}

std::ostream& operator<<(std::ostream& os, const Mesh &m)
{
	for (size_t i = 0; i < m._elements.size(); i++) {
		os << *(m._elements[i]) << "\n";
	}
	return os;
}

}


