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
		eslocal *eSubPartition = getPartition(_partPtrs[i], _partPtrs[i + 1], fixPoints);

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

void Mesh::computeCorners(eslocal number, bool vertex, bool edges, bool faces, bool averageEdges, bool averageFaces)
{
	if (parts() < 1) {
		std::cerr << "Internal error: _partPtrs.size()\n";
		exit(EXIT_FAILURE);
	}
	if (parts() == 1 || (!vertex && !edges && !faces && !averageEdges && !averageFaces)) {
		return;
	}

	// node to element
	std::vector<std::vector<eslocal> > nodesElements(_coordinates.clusterSize());

	for (size_t i = 0; i < parts(); i++) {
		for (eslocal j = _partPtrs[i]; j < _partPtrs[i + 1]; j++) {
			for (size_t k = 0; k < _elements[j]->size(); k++) {
				// add node's adjacent element
				nodesElements[_coordinates.clusterIndex(_elements[j]->node(k), i)].push_back(j);
			}
		}
	}

	// vector of faces
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
				if (findSubdomains(nodesElements, face, _partPtrs, 1, eSub) == j) {
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
	Mesh cfm;
	Element *el;
	eslocal index = 0;
	for (size_t i = 0; i < _coordinates.clusterSize(); i++) {
		if (projection[i] == 1) {
			cfm.coordinates().add(_coordinates[i], index, i);
			projection[i] = index++;
		}
	}
	for (size_t i = 0; i < commonFaces.size(); i++) {
		for (size_t j = 0; j < commonFaces[i].size(); j++) {
			commonFaces[i][j] = projection[commonFaces[i][j]];
		}
	}
	cfm._elements.reserve(commonFaces.size());

	// create mesh
	cfm._partPtrs.clear();
	cfm._partPtrs.push_back(0);
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
					cfm._elements.push_back(el);
				}
			}
			if (cfm._elements.size() > cfm._partPtrs.back()) {
				cfm._partPtrs.push_back(cfm._elements.size());
			}
		}
	}
	cfm.remapElementsToSubdomain();

	if (faces) {
		cfm.computeFixPoints(number);
		for (size_t p = 0; p < cfm.parts(); p++) {
			for (size_t i = 0; i < cfm.getFixPoints()[p].size(); i++) {
				_subdomainBoundaries.setCorner(cfm.coordinates().globalIndex(cfm.getFixPoints()[p][i], p));
			}
		}
	}

	if (!edges && !vertex && !averageEdges && !averageFaces) {
		return;
	}

	Mesh clm;
	std::vector<std::vector<std::vector<eslocal> > > nodesFaces(cfm.parts());
	for (size_t p = 0; p < cfm.parts(); p++) {
		nodesFaces[p].resize(cfm._coordinates.localSize(p));
		for (eslocal e = cfm._partPtrs[p]; e < cfm._partPtrs[p + 1]; e++) {
			for (size_t n = 0; n < cfm._elements[e]->size(); n++) {
				// add node's adjacent element
				nodesFaces[p][cfm._elements[e]->node(n)].push_back(e);
			}
		}
	}

	std::set<eslocal> borderLine;

	std::vector<std::pair<eslocal, eslocal> > outerFaces;
	std::vector<char> nSubdomains(cfm.coordinates().clusterSize() * cfm.parts(), 0);
	std::vector<eslocal> points(cfm.coordinates().clusterSize(), 0);
	clm._partPtrs.clear();
	clm._partPtrs.push_back(0);
	for (size_t p = 0; p < cfm.parts(); p++) {
		for (size_t e = cfm._partPtrs[p]; e < cfm._partPtrs[p + 1]; e++) {
			for (size_t f = 0; f < cfm._elements[e]->faces(); f++) {
				std::vector<eslocal> face = cfm._elements[e]->getFace(f);
				if (isOuterFace(nodesFaces[p], face)) {
					if (averageFaces) {
						borderLine.insert(cfm.coordinates().clusterIndex(face[0], p));
						borderLine.insert(cfm.coordinates().clusterIndex(face[1], p));
					}
					face[0] = cfm.coordinates().clusterIndex(face[0], p);
					face[1] = cfm.coordinates().clusterIndex(face[1], p);
					if (face[0] < face[1]) {
						outerFaces.push_back(std::pair<eslocal, eslocal>(face[0], face[1]));
					} else {
						outerFaces.push_back(std::pair<eslocal, eslocal>(face[1], face[0]));
					}
					nSubdomains[face[0] * cfm.parts() + p] = 1;
					nSubdomains[face[1] * cfm.parts() + p] = 1;
					points[face[0]] = 1;
					points[face[1]] = 1;
				}
			}
		}
	}
	std::sort(outerFaces.begin(), outerFaces.end());
	auto it = std::unique(outerFaces.begin(), outerFaces.end());
	outerFaces.resize(std::distance(outerFaces.begin(), it));

	// faces averaging
	if (averageFaces) {
		if (!cfm.getFixPoints().size() || !cfm.getFixPoints()[0].size()) {
			cfm.computeFixPoints(1);
		}
		for (size_t p = 0; p < cfm.parts(); p++) {
			eslocal corner = cfm.coordinates().globalIndex(cfm.getFixPoints()[p][0], p);
			_subdomainBoundaries.setCorner(corner);
			std::set<eslocal> aPoints;
			for (size_t e = cfm._partPtrs[p]; e < cfm._partPtrs[p + 1]; e++) {
				for (size_t n = 0; n < cfm._elements[e]->size(); n++) {
					eslocal node = cfm.coordinates().clusterIndex(cfm._elements[e]->node(n), p);
					if (borderLine.find(node) == borderLine.end()) {
						aPoints.insert(cfm.coordinates().globalIndex(node));
					}
				}
			}
			std::vector<eslocal>& averaging = _subdomainBoundaries.averaging(corner);
			averaging.insert(averaging.begin(), aPoints.begin(), aPoints.end());
		}
	}

	eslocal linePoints = 0;
	for (size_t i = 0; i < cfm._coordinates.clusterSize(); i++) {
		if (points[i] == 1) {
			clm.coordinates().add(cfm._coordinates[i], linePoints, i);
			points[i] = linePoints++;
		}
	}

	std::vector<bool> mask(cfm.coordinates().clusterSize(), false);
	size_t maskCounter = 0;
	for (size_t i = 0; i < mask.size(); i++) {
		char *s = nSubdomains.data() + i * cfm.parts();
		char *e = nSubdomains.data() + (i + 1) * cfm.parts();
		if (std::all_of(s, e, [](char flag) { return flag == 0; })) {
			mask[i] = true;
			maskCounter++;
		}
	}

	size_t begin = 0;
	std::vector<eslocal> tmpPair(2);
	std::vector<bool> sausageEdges;
	while (maskCounter < cfm.coordinates().clusterSize()) {
		while (mask[begin]) {
			begin++;
		}
		mask[begin] = true;
		char *s = nSubdomains.data() + begin * cfm.parts();
		char *e = nSubdomains.data() + (begin + 1) * cfm.parts();
		for (size_t i = begin; i < mask.size(); i++) {
			char *t = nSubdomains.data() + i * cfm.parts();
			if (std::equal(s, e, t)) {
				mask[i] = true;
				maskCounter++;
			}
		}
		bool innerLine = std::count(s, e, 1) > 1;
		std::set<eslocal> vertices;
		bool flag = true;
		for (size_t i = 0; i < outerFaces.size(); i++) {
			char *p0s = nSubdomains.data() + outerFaces[i].first * cfm.parts();
			char *p1s = nSubdomains.data() + outerFaces[i].second * cfm.parts();
			char *p0e = nSubdomains.data() + (outerFaces[i].first + 1) * cfm.parts();
			char *p1e = nSubdomains.data() + (outerFaces[i].second + 1) * cfm.parts();
			if (std::equal(s, e, p0s) && std::equal(s, e, p1s)) {
				tmpPair[0] = clm.coordinates().clusterIndex(outerFaces[i].first);
				tmpPair[1] = clm.coordinates().clusterIndex(outerFaces[i].second);
				if (averageEdges) {
					esglobal p1 = cfm.coordinates().globalIndex(outerFaces[i].first);
					esglobal p2 = cfm.coordinates().globalIndex(outerFaces[i].second);
					auto &dx = _coordinates.property(DIRICHLET_X).values();
					auto &dy = _coordinates.property(DIRICHLET_Y).values();
					auto &dz = _coordinates.property(DIRICHLET_Z).values();
					auto &cb = _clusterBoundaries;
					if (dx.find(p1) == dx.end() && dx.find(p2) == dx.end()
							&& dy.find(p1) == dy.end() && dy.find(p2) == dy.end()
							&& dz.find(p1) == dz.end() && dz.find(p2) == dz.end()
							&& cb[p1].size() == 1 && cb[p2].size() == 1) {

						clm._elements.push_back(new Line(tmpPair.data()));
					}
				} else {
					clm._elements.push_back(new Line(tmpPair.data()));
				}
			}
			if (std::equal(s, e, p0s) != std::equal(s, e, p1s)) {
				if (std::count(p0s, p0e, 1) > std::count(p1s, p1e, 1)) {
					vertices.insert(clm.coordinates().clusterIndex(outerFaces[i].first));
				} else {
					vertices.insert(clm.coordinates().clusterIndex(outerFaces[i].second));
				}
			}
		}
		for (eslocal e = clm._elements.size() - 1; e >= clm._partPtrs.back(); e--) {
			if (vertices.find(clm._elements[e]->node(0)) != vertices.end()
					|| vertices.find(clm._elements[e]->node(1)) != vertices.end()) {

				// erase elements with vertex corners
				clm._elements.erase(clm._elements.begin() + e);
			}
		}
		if (vertices.size() > 2) {
			// non-continuous edge -> separate it
			std::vector<std::pair<eslocal, eslocal> > neigh(clm.coordinates().clusterSize(), std::pair<eslocal, eslocal>(-1, -1));
			std::vector<int> color(clm._elements.size() - clm._partPtrs.back());
			auto next = [&] (eslocal node) -> eslocal {
				if (neigh[node].first != -1 && color[neigh[node].first - clm._partPtrs.back()] == 0) {
					return neigh[node].first;
				};
				if (neigh[node].second != -1 && color[neigh[node].second - clm._partPtrs.back()] == 0) {
					return neigh[node].second;
				};
				return -1;
			};
			std::function<void(eslocal, int)> traverse = [&] (eslocal element, int lColor) {
				color[element - clm._partPtrs.back()] = lColor;
				eslocal e = next(clm._elements[element]->node(0));
				if (e != -1 && color[e - clm._partPtrs.back()] == 0) {
					traverse(e, lColor);
				}
				e = next(clm._elements[element]->node(1));
				if (e != -1 && color[e - clm._partPtrs.back()] == 0) {
					traverse(e, lColor);
				}
			};

			for (eslocal e = clm._partPtrs.back(); e < clm._elements.size(); e++) {
				if (neigh[clm._elements[e]->node(0)].first == -1) {
					neigh[clm._elements[e]->node(0)].first = e;
				} else {
					neigh[clm._elements[e]->node(0)].second = e;
				}
				if (neigh[clm._elements[e]->node(1)].first == -1) {
					neigh[clm._elements[e]->node(1)].first = e;
				} else {
					neigh[clm._elements[e]->node(1)].second = e;
				}
			}

			int colors = 0;
			for (eslocal e = clm._partPtrs.back(), i = 0; e < clm._elements.size(); e++, i++) {
				if (color[i] == 0) {
					traverse(e, ++colors);
				}
			}
			std::vector<Element*> tmp(clm._elements.begin() + clm._partPtrs.back(), clm._elements.end());
			clm._elements.resize(clm._partPtrs.back());
			for (int c = 1; c <= colors; c++) {
				for (size_t i = 0; i < color.size(); i++) {
					if (color[i] == c) {
						clm._elements.push_back(tmp[i]);
					}
				}
				clm._partPtrs.push_back(clm._elements.size());
			}
		}

		if (clm._partPtrs.back() < clm._elements.size()) {
			clm._partPtrs.push_back(clm._elements.size());
		}

		if (vertex) {
			sausageEdges.push_back(!vertices.size());
			for (auto it = vertices.begin(); it != vertices.end(); ++it) {
				eslocal fPoint = clm.coordinates().globalIndex(*it);
				_subdomainBoundaries.setCorner(cfm.coordinates().globalIndex(fPoint));
			}
		}
	}

	clm.remapElementsToSubdomain();
	for (size_t i = 0; i < sausageEdges.size(); i++) {
		if (sausageEdges[i] && 0) { 
			size_t SAUSAGE_CORNERS = 4;
			eslocal *eSubPartition = clm.getPartition(clm._partPtrs[i], clm._partPtrs[i + 1], SAUSAGE_CORNERS);

			for (eslocal j = 0; j < SAUSAGE_CORNERS; j++) {
				eslocal sCorner = clm.getCentralNode(clm._partPtrs[i], clm._partPtrs[i + 1], eSubPartition, i, j);
				sCorner = clm.coordinates().globalIndex(sCorner, i);
				_subdomainBoundaries.setCorner(cfm.coordinates().globalIndex(sCorner));
			}

			delete[] eSubPartition;
		}
	}

	if (averageEdges) {
		clm.computeFixPoints(1);
		for (size_t p = 0; p < clm.parts(); p++) {
			for (size_t i = 0; i < clm.getFixPoints()[p].size(); i++) {
				eslocal corner = clm.coordinates().globalIndex(clm.getFixPoints()[p][i], p);
				corner = cfm.coordinates().globalIndex(corner);
				_subdomainBoundaries.setCorner(corner);
				std::set<eslocal> aPoints;
				for (size_t e = clm._partPtrs[p]; e < clm._partPtrs[p + 1]; e++) {
					eslocal p0 = clm.coordinates().globalIndex(clm._elements[e]->node(0), p);
					eslocal p1 = clm.coordinates().globalIndex(clm._elements[e]->node(1), p);
					p0 = cfm.coordinates().globalIndex(p0);
					p1 = cfm.coordinates().globalIndex(p1);
					if (!_subdomainBoundaries.isCorner(p0)) {
						aPoints.insert(p0);
					}
					if (!_subdomainBoundaries.isCorner(p1)) {
						aPoints.insert(p1);
					}
				}
				std::vector<eslocal>& averaging = _subdomainBoundaries.averaging(corner);
				averaging.insert(averaging.begin(), aPoints.begin(), aPoints.end());
			}
		}
		return;
	}
	if (edges) {
		clm.computeFixPoints(number);
		for (size_t p = 0; p < clm.parts(); p++) {
			for (size_t i = 0; i < clm.getFixPoints()[p].size(); i++) {
				eslocal fPoint = clm.coordinates().globalIndex(clm.getFixPoints()[p][i], p);
				_subdomainBoundaries.setCorner(cfm.coordinates().globalIndex(fPoint));
			}
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


