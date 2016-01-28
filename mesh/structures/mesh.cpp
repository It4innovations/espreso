#include "mesh.h"

#include "esinput.h"
#include "esoutput.h"

using namespace mesh;

Mesh::Mesh():_elements(0), _fixPoints(0)
{
	_partPtrs.resize(2);
	_partPtrs[0] = 0;
	_partPtrs[1] = 0;
}

void Mesh::partitiate(size_t parts)
{
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
				_subdomainBoundaries[_coordinates.clusterIndex(_elements[e]->node(n), p)].insert(p);
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
		eslog::error << "An input for METIS procedure is incorrect.\n";
		exit(EXIT_FAILURE);
	case METIS_ERROR_MEMORY:
		eslog::error << "There is not enough memory for compute a partition.\n";
		exit(EXIT_FAILURE);
	case METIS_ERROR:
		eslog::error << "METIS fail computation.\n";
		exit(EXIT_FAILURE);
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

static eslocal findSubdomains(
		std::vector<std::vector<eslocal> > &nodesElements,
		std::vector<eslocal> &face,
		const std::vector<eslocal> &partPtrs,
		eslocal differentParts,
		std::vector<eslocal> &subdomains)
{
	subdomains.clear();
	eslocal NOT_ON_BOUNDARY = -1;
	std::vector<eslocal> result(nodesElements[face[0]]);
	std::vector<eslocal>::iterator it = result.end();

	for (size_t i = 1; i < face.size(); i++) {
		std::vector<eslocal> tmp(result.begin(), it);
		it = std::set_intersection(tmp.begin(), tmp.end(),
				nodesElements[face[i]].begin(), nodesElements[face[i]].end(),
				result.begin());
		if (it - result.begin() == 1) {
			return NOT_ON_BOUNDARY;
		}
	}

	size_t r = 0;
	for (size_t p = 1; p < partPtrs.size(); p++) {
		for ( ; r < it - result.begin() && result[r] < partPtrs[p]; r++) {
			if (!subdomains.size() || subdomains.back() != p - 1) {
				subdomains.push_back(p - 1);
			}
		}
	}

	return (subdomains.size() > differentParts) ? result[0] : NOT_ON_BOUNDARY;
}

void Mesh::getSurface(Mesh &surface) const
{
	// vector of faces in all parts
	std::vector<std::vector<std::vector<eslocal> > > faces(parts());
	// number of elements in all parts
	std::vector<size_t> elementsCount(parts(), 0);

	if (parts() < 1) {
		std::cerr << "Internal error: _partPtrs.size()\n";
		exit(EXIT_FAILURE);
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

	// create surface mesh
	surface._partPtrs.push_back(0); //(surface._elements.size());
	for (size_t i = 0; i + 1 < _partPtrs.size(); i++) {
		for (size_t j = 0; j < faces[i].size(); j++) {
			std::vector<eslocal> &face = faces[i][j];
			if (face.size() == 3) {
				surface._elements.push_back(new Triangle(&face[0]));
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
					surface._elements.push_back(new Triangle(&face[0]));
					face[1] = face[0];
					surface._elements.push_back(new Triangle(&face[1]));
				} else {
					surface._elements.push_back(new Triangle(&face[1]));
					face[2] = face[3];
					surface._elements.push_back(new Triangle(&face[0]));
				}
			}
		}
		surface._partPtrs.push_back(surface._elements.size());
	}

	surface.remapElementsToSubdomain();
	surface.computeFixPoints(0);
	surface.computeBoundaries();
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

static std::vector<std::vector<eslocal> > getNodeToElementsMap(const Mesh &mesh)
{
	auto &partPtrs = mesh.getPartition();
	auto &elements = mesh.getElements();

	std::vector<std::vector<eslocal> > map(mesh.coordinates().clusterSize());

	for (size_t i = 0; i < mesh.parts(); i++) {
		for (eslocal j = partPtrs[i]; j < partPtrs[i + 1]; j++) {
			for (size_t k = 0; k < elements[j]->size(); k++) {
				// add node's adjacent element
				map[mesh.coordinates().clusterIndex(elements[j]->node(k), i)].push_back(j);
			}
		}
	}

	return map;
}

void Mesh::computeCommonFaces(Mesh &mesh)
{
	mesh._coordinates.clear();
	mesh._elements.clear();
	mesh._partPtrs.clear();

	std::vector<std::vector<eslocal> > nodeToElements = getNodeToElementsMap(*this);
	std::vector<std::vector<eslocal> > commonFaces;
	std::vector<eslocal> eSub;
	std::vector<bool> subdomains;
	std::vector<eslocal> projection(_coordinates.clusterSize(), 0);

	for (eslocal i = 0; i < parts(); i++) {
		// compute number of elements and fill used nodes
		for (eslocal j = _partPtrs[i]; j < _partPtrs[i + 1]; j++) {
			for (size_t k = 0; k < _elements[j]->faces(); k++) {
				std::vector<eslocal> face = _elements[j]->getFace(k);
				for (size_t f = 0; f < face.size(); f++) {
					face[f] = _coordinates.clusterIndex(face[f], i);
				}
				if (findSubdomains(nodeToElements, face, _partPtrs, 1, eSub) == j) {
					subdomains.resize(subdomains.size() + parts(), false);
					subdomains[subdomains.size() - parts() + eSub[0]] = true;
					subdomains[subdomains.size() - parts() + eSub[1]] = true;
					commonFaces.push_back(face);
					for (size_t f = 0; f < face.size(); f++) {
						projection[face[f]] = 1;
					}
				}
			}
		}
	}

	// create mesh from common faces
	Element *el;
	eslocal index = 0;
	for (size_t i = 0; i < _coordinates.clusterSize(); i++) {
		if (projection[i] == 1) {
			mesh.coordinates().add(_coordinates[i], index, i);
			projection[i] = index++;
		}
	}
	for (size_t i = 0; i < commonFaces.size(); i++) {
		for (size_t j = 0; j < commonFaces[i].size(); j++) {
			commonFaces[i][j] = projection[commonFaces[i][j]];
		}
	}
	mesh._elements.reserve(commonFaces.size());

	// create mesh
	mesh._partPtrs.clear();
	mesh._partPtrs.push_back(0);
	for (size_t i = 0; i < parts(); i++) {
		for (size_t j = i + 1; j < parts(); j++) {
			for (size_t e = 0; e < commonFaces.size(); e++) {
				if (subdomains[e * parts() + i] && subdomains[e * parts() + j]) {
					if (commonFaces[e].size() == 3) {
						el = new Triangle(commonFaces[e].data());
					}
					if (commonFaces[e].size() == 4) {
						el = new Square(commonFaces[e].data());
					}
					mesh._elements.push_back(el);
				}
			}
			if (mesh._elements.size() > mesh._partPtrs.back()) {
				mesh._partPtrs.push_back(mesh._elements.size());
				mesh.makePartContinuous(mesh.parts() - 1);
			}
		}
	}
	mesh.remapElementsToSubdomain();
	mesh.computeFixPoints(0);
}

void Mesh::computeBorderLinesAndVertices(const Mesh &faces,std::vector<char> &border, Mesh &lines, std::set<eslocal> &vertices)
{
	lines._elements.clear();
	lines._partPtrs.clear();
	lines._coordinates.clear();
	border.clear();
	border.resize(faces.coordinates().clusterSize(), 0);

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
						border[line[n]] = 1;
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
						std::cerr << "MESH ERROR: unknown line type.\n";
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

	// check whether a point has dirichler condition
	auto has_dirichlet = [&] (eslocal i) {
		eslocal index = faces.coordinates().globalIndex(i);
		auto &dx = coordinates().property(DIRICHLET_X).values();
		auto &dy = coordinates().property(DIRICHLET_Y).values();
		auto &dz = coordinates().property(DIRICHLET_Z).values();
		return dx.find(index) != dx.end() || dy.find(index) != dy.end() || dz.find(index) != dz.end();
	};

	// check whether a point is on cluster boundary
	auto on_cluster_boundary = [&] (eslocal i) {
		eslocal index = faces.coordinates().globalIndex(i);
		auto &cb = clusterBoundaries();
		return cb[index].size() > 1;
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
		for (size_t i = 0; i < commonLines.size(); i++) {
			eslocal start = std::get<0>(commonLines[i]);
			eslocal mid = std::get<1>(commonLines[i]);
			eslocal end = (std::get<2>(commonLines[i]) == -1) ? mid : std::get<2>(commonLines[i]);

			if (same_subdomains(begin, start) && same_subdomains(begin, end)) {
				if (!has_dirichlet(start) && !has_dirichlet(end) && !on_cluster_boundary(start) && !on_cluster_boundary(end)) {
					eslocal tmp[3] = { points[start], points[mid], points[end] };
					if (std::get<2>(commonLines[i]) == -1) {
						lines._elements.push_back(new Line(tmp));
					} else {
						lines._elements.push_back(new Line(tmp));
					}
				} else {
					if (has_dirichlet(start) || on_cluster_boundary(end)) {
						vertices.insert(points[start]);
					} else {
						vertices.insert(points[end]);
					}
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

void Mesh::checkCycle(size_t part, std::set<eslocal> &vertices)
{
	std::vector<int> counter(_coordinates.localSize(part), 0);
	for (size_t e = _partPtrs[part]; e < _partPtrs[part + 1]; e++) {
		counter[_elements[e]->node(0)]++;
		counter[_elements[e]->node(_elements[e]->size() - 1)]++;
	}
	if (std::all_of(counter.begin(), counter.end(), [] (int count) { return count % 2 == 0; })) {
		int max = (_partPtrs[part + 1] - _partPtrs[part]) / 5 + 1;
		// TODO: try this:
		//int max = _partPtrs[p + 1] - _partPtrs[p];
		size_t corners = std::min(max, 4);
		eslocal *eSubPartition = getPartition(_partPtrs[part], _partPtrs[part + 1], corners);
		for (eslocal j = 0; j < corners; j++) {
			eslocal sCorner = getCentralNode(_partPtrs[part], _partPtrs[part + 1], eSubPartition, part, j);
			vertices.insert(_coordinates.clusterIndex(sCorner, part));
		}
		delete[] eSubPartition;
	}
}


void Mesh::computeCorners(eslocal number, bool vertices, bool edges, bool faces, bool averageEdges, bool averageFaces)
{
	if (parts() < 1) {
		std::cerr << "Internal error: _partPtrs.size()\n";
		exit(EXIT_FAILURE);
	}
	if (parts() == 1 || (!vertices && !edges && !faces && !averageEdges && !averageFaces)) {
		return;
	}

	Mesh commonFaces;
	Mesh commonLines;
	std::set<eslocal> commonVertices;
	std::vector<char> commonFacesBorner;

	computeCommonFaces(commonFaces);
	computeBorderLinesAndVertices(commonFaces, commonFacesBorner, commonLines, commonVertices);


	auto faceToCluster = [&] (eslocal index, eslocal part) {
		return commonFaces.coordinates().globalIndex(index, part);
	};

	auto lineToCluster = [&] (eslocal index, eslocal part) {
		index = commonLines.coordinates().globalIndex(index, part);
		return commonFaces.coordinates().globalIndex(index);
	};

	if (vertices) {
		for (auto it = commonVertices.begin(); it != commonVertices.end(); ++it) {
			_subdomainBoundaries.setCorner(commonFaces.coordinates().globalIndex(*it));
		}
		if (!edges && !averageEdges) {
			std::set<eslocal> cycleVertices;
			for (size_t p = 0; p < commonLines.parts(); p++) {
				commonLines.checkCycle(p, cycleVertices);
			}
			for (auto it = cycleVertices.begin(); it != cycleVertices.end(); ++it) {
				_subdomainBoundaries.setCorner(commonFaces.coordinates().globalIndex(*it));
			}
		}
	}

	if (edges && !averageEdges) {
		commonLines.computeFixPoints(number);
		for (size_t p = 0; p < commonLines.parts(); p++) {
			auto &fixPoints = commonLines.getFixPoints();
			for (size_t i = 0; i < fixPoints[p].size(); i++) {
				_subdomainBoundaries.setCorner(lineToCluster(fixPoints[p][i], p));
			}
		}
	}

	if (averageEdges) {
		commonLines.computeFixPoints(1);
		for (size_t p = 0; p < commonLines.parts(); p++) {
			eslocal corner = lineToCluster(commonLines._fixPoints[p][0], p);
			_subdomainBoundaries.setCorner(corner);

			std::set<eslocal> aPoints;
			for (size_t e = commonLines._partPtrs[p]; e < commonLines._partPtrs[p + 1]; e++) {
				for (size_t n = 0; n < commonLines._elements[e]->size(); n++) {
					eslocal point = lineToCluster(commonLines._elements[e]->node(n), p);
					if (!_subdomainBoundaries.isCorner(point)) {
						aPoints.insert(point);
					}
				}
			}

			std::vector<eslocal>& averaging = _subdomainBoundaries.averaging(corner);
			averaging.insert(averaging.begin(), aPoints.begin(), aPoints.end());
		}
	}

	size_t faceCorners = 0;
	if (faces) {
		faceCorners += number;
	}
	if (averageFaces) {
		faceCorners++;
	}
	commonFaces.computeFixPoints(faceCorners);

	if (faces) {
		for (size_t p = 0; p < commonFaces.parts(); p++) {
			auto &fixPoints = commonFaces.getFixPoints();
			for (size_t i = 0; i < fixPoints[p].size(); i++) {
				_subdomainBoundaries.setCorner(faceToCluster(fixPoints[p][i], p));
			}
		}
	}
	if (averageFaces) {
		for (size_t p = 0; p < commonFaces.parts(); p++) {
			eslocal corner = -1;
			for (size_t i = 0; i < commonFaces._fixPoints[p].size(); i++) {
				eslocal point = commonFaces.coordinates().clusterIndex(commonFaces._fixPoints[p][i], p);
				if (commonFacesBorner[point] == 0) {
					corner = commonFaces.coordinates().globalIndex(point);
					break;
				}
			}
			if (corner == -1) {
				for (size_t i = 0; i < commonFaces.coordinates().localSize(p); i++) {
					eslocal point = commonFaces.coordinates().clusterIndex(i, p);
					if (commonFacesBorner[point] == 0) {
						corner = commonFaces.coordinates().globalIndex(point);
						break;
					}
				}
			}

			if (corner == -1) {
				std::cout << "AVERAGE FACES WARNING: There is no point inside common faces.\n";
				break;
			}

			_subdomainBoundaries.setCorner(corner);

			std::set<eslocal> aPoints;
			for (size_t e = commonFaces._partPtrs[p]; e < commonFaces._partPtrs[p + 1]; e++) {
				for (size_t n = 0; n < commonFaces._elements[e]->size(); n++) {
					eslocal gPoint = faceToCluster(commonFaces._elements[e]->node(n), p);
					eslocal cPoint = commonFaces.coordinates().clusterIndex(commonFaces._elements[e]->node(n), p);
					if (!_subdomainBoundaries.isCorner(gPoint) && !commonFacesBorner[cPoint]) {
						aPoints.insert(gPoint);
					}
				}
			}

			std::vector<eslocal>& averaging = _subdomainBoundaries.averaging(corner);
			averaging.insert(averaging.begin(), aPoints.begin(), aPoints.end());
		}
	}
}

void Mesh::remapElementsToSubdomain()
{
	std::vector<eslocal> nodeMap(_coordinates.clusterSize(), -1);
	_coordinates._clusterIndex.clear();
	_coordinates._clusterIndex.resize(parts());

	for (size_t p = 0; p < parts(); p++) {
		std::fill(nodeMap.begin(), nodeMap.end(), -1);

		// Compute mask of nodes
		for (eslocal e = _partPtrs[p]; e < _partPtrs[p + 1]; e++) {
			for (size_t n = 0; n < _elements[e]->size(); n++) {
				nodeMap[_elements[e]->node(n)] = 1;
			}
		}

		// re-index nodes
		eslocal nSize = 0;
		for (size_t c = 0; c < _coordinates.clusterSize(); c++) {
			if (nodeMap[c] == 1) {
				nodeMap[c] = nSize++;
			}
		}

		for (eslocal e = _partPtrs[p]; e < _partPtrs[p + 1]; e++) {
			for (eslocal n = 0; n < _elements[e]->size(); n++) {
				_elements[e]->node(n) = nodeMap[_elements[e]->node(n)];
			}
		}

		_coordinates._clusterIndex[p].reserve(nSize);
		for (size_t c = 0; c < _coordinates.clusterSize(); c++) {
			if (nodeMap[c] >= 0) {
				_coordinates._clusterIndex[p].push_back(c);
			}
		}
	}
}

void Mesh::remapElementsToCluster()
{
	for (eslocal p = 0; p < this->parts(); p++) {
		for (eslocal e = _partPtrs[p]; e < _partPtrs[p + 1]; e++) {
			for (eslocal n = 0; n < _elements[e]->size(); n++) {
				_elements[e]->node(n) = _coordinates.clusterIndex(_elements[e]->node(n), p);
			}
		}
	}
}

std::ostream& mesh::operator<<(std::ostream& os, const Mesh &m)
{
	for (size_t i = 0; i < m._elements.size(); i++) {
		os << *(m._elements[i]) << "\n";
	}
	return os;
}


