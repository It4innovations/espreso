#include "mesh.h"

#include "esinput.h"
#include "esoutput.h"

using namespace mesh;

Mesh::Mesh(int rank, int size):_elements(0), _fixPoints(0), _rank(rank), _size(size)
{
	_partPtrs.resize(2);
	_partPtrs[0] = 0;
	_partPtrs[1] = 0;
}

void Mesh::partitiate(eslocal parts, eslocal fixPoints)
{
	if (this->parts()) {
		// reset elements node indices
		for (eslocal p = 0; p < this->parts(); p++) {
			for (eslocal e = _partPtrs[p]; e < _partPtrs[p + 1]; e++) {
				for (eslocal n = 0; n < _elements[e]->size(); n++) {
					_elements[e]->node(n) = _coordinates.clusterIndex(_elements[e]->node(n), p);
				}
			}
		}
	}
	_partPtrs.resize(parts + 1);
	_coordinates.localClear();
	_coordinates.localResize(parts);

	// Call METIS to get partition of a whole mesh
	eslocal *elementPartition = getPartition(0, _elements.size(), parts);

	// Rearrange mesh's elements
	partitiate(elementPartition);

	delete[] elementPartition;

	if (fixPoints > 0) {
		computeFixPoints(fixPoints);
	} else {
		_fixPoints.resize(parts, std::vector<eslocal>());
	}

	computeBoundaries();
}

void Mesh::computeBoundaries()
{
	_subdomainBoundaries.clear();
	_subdomainBoundaries.resize(_coordinates.size());

	for (size_t p = 0; p < parts(); p++) {
		for (eslocal e = _partPtrs[p]; e < _partPtrs[p + 1]; e++) {
			for (size_t n = 0; n < _elements[e]->size(); n++) {
				_subdomainBoundaries[_coordinates.clusterIndex(_elements[e]->node(n), p)].insert(p);
			}
		}
	}
}

void Mesh::computeFixPoints(eslocal fixPoints)
{
	_fixPoints.clear();
	_fixPoints.resize(parts(), std::vector<eslocal>(fixPoints));

#ifndef DEBUG
	cilk_for (eslocal i = 0; i < parts(); i++) {
#else
	for (eslocal i = 0; i < parts(); i++) {
#endif
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
		e[index + 1] = e[index] + _elements[i]->size();
	}

	// create array of nodes
	n = new eslocal[e[eSize]];
	// number of common nodes to be neighbor
	ncommon = 4;
	for (eslocal i = first, index = 0; i < last; i++, index++) {
		_elements[i]->fillNodes(n + e[index]);
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

void Mesh::partitiate(eslocal *ePartition)
{
	_partPtrs[0] = 0;

	Element *e;
	eslocal p;
	for (size_t part = 0; part < _partPtrs.size() - 1; part++) {
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
		computeLocalIndices(part);
	}
}

void Mesh::computeLocalIndices(size_t part)
{
	std::vector<eslocal> nodeMap(_coordinates.clusterSize(), -1);

	// Compute mask of nodes
	for (eslocal e = _partPtrs[part]; e < _partPtrs[part + 1]; e++) {
		for (size_t n = 0; n < _elements[e]->size(); n++) {
			nodeMap[_elements[e]->node(n)] = 1;
		}
	}

	// re-index nodes
	eslocal nSize = 0;
	for (eslocal k = 0; k < _coordinates.clusterSize(); k++) {
		if (nodeMap[k] == 1) {
			nodeMap[k] = nSize++;
		}
	}

	for (eslocal e = _partPtrs[part]; e < _partPtrs[part + 1]; e++) {
		for (eslocal n = 0; n < _elements[e]->size(); n++) {
			_elements[e]->node(n) = nodeMap[_elements[e]->node(n)];
		}
	}

	_coordinates.computeLocal(part, nodeMap, nSize);
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

void Mesh::checkMETISResult(eslocal result) const
{
	switch (result) {
	case METIS_ERROR_INPUT:
		fprintf(stderr, "An input for METIS procedure is incorrect.\n");
		exit(EXIT_FAILURE);
	case METIS_ERROR_MEMORY:
		fprintf(stderr, "There is not enough memory for compute a partition.\n");
		exit(EXIT_FAILURE);
	case METIS_ERROR:
		fprintf(stderr, "METIS fail computation.\n");
		exit(EXIT_FAILURE);
	}
}

void Mesh::checkMKLResult(eslocal result) const
{
	switch (result) {
	case 0:
		return;
	default:
		std::cerr << "MKL error: " << result << ".\n";
		exit(EXIT_FAILURE);
	}
}

Mesh::~Mesh()
{
	for (size_t i = 0; i < _elements.size(); i++) {
		delete _elements[i];
	}
}

void Mesh::saveNodeArray(eslocal *nodeArray, size_t part)
{
	size_t p = 0;
	for (eslocal i = _partPtrs[part]; i < _partPtrs[part + 1]; i++) {
		_elements[i]->fillNodes(nodeArray + p);
		p += _elements[i]->size();
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

static eslocal isOnBoundary(
		std::vector<std::vector<eslocal> > &nodesElements,
		std::vector<eslocal> &face,
		const std::vector<eslocal> &partPtrs,
		eslocal differentParts)
{
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

	eslocal counter = 0;
	eslocal minPart = result[0];
	for (size_t r = 1; r < it - result.begin(); r++) {
		for (size_t i = 1; i < partPtrs.size() - 1; i++) {
			if (result[r - 1] < partPtrs[i] && partPtrs[i] <= result[r]) {
				if (minPart > result[r]) {
					minPart = result[r];
				}
				counter++;
				break;
			}
		}
	}
	return (counter >= differentParts) ? minPart : NOT_ON_BOUNDARY;
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

void Mesh::getSurface(SurfaceMesh &surface) const
{
	// vector of faces in all parts
	std::vector<std::vector<std::vector<eslocal> > > faces(parts());
	// number of elements in all parts
	std::vector<size_t> elementsCount(parts(), 0);

	if (parts() < 1) {
		std::cerr << "Internal error: _partPtrs.size()\n";
		exit(EXIT_FAILURE);
	}
#ifndef DEBUG
	cilk_for (size_t i = 0; i < parts(); i++) {
#else
	for (size_t i = 0; i < parts(); i++) {
#endif
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
		surface.computeLocalIndices(surface._partPtrs.size() - 2);

	}

	surface.computeBoundaries();

}

void Mesh::getCommonFaces(CommonFacesMesh &commonFaces) const
{
	// vector of faces in all parts
	std::vector<std::vector<std::vector<eslocal> > > faces(parts());
	// number of elements in all parts
	std::vector<size_t> elementsCount(parts(), 0);

	if (parts() < 1) {
		std::cerr << "Internal error: _partPtrs.size()\n";
		exit(EXIT_FAILURE);
	}

	std::vector<std::vector<eslocal> > nodesElements(_coordinates.size());
	for (size_t i = 0; i < parts(); i++) {
		// Compute nodes' adjacent elements
		for (eslocal j = _partPtrs[i]; j < _partPtrs[i + 1]; j++) {
			for (size_t k = 0; k < _elements[j]->size(); k++) {
				nodesElements[_coordinates.clusterIndex(_elements[j]->node(k), i)].push_back(j);
			}
		}
	}

#ifndef DEBUG
	cilk_for (size_t i = 0; i < parts(); i++) {
#else
	for (size_t i = 0; i < parts(); i++) {
#endif
		// compute number of elements and fill used nodes
		for (eslocal j = _partPtrs[i]; j < _partPtrs[i + 1]; j++) {
			for (size_t k = 0; k < _elements[j]->faces(); k++) {
				std::vector<eslocal> face = _elements[j]->getFace(k);
				for (size_t f = 0; f < face.size(); f++) {
					face[f] = _coordinates.clusterIndex(face[f], i);
				}
				if (isOnBoundary(nodesElements, face, _partPtrs, 1) >= 0) {
					faces[i].push_back(face);
					elementsCount[i] += 1;
				}
			}
		}
	}

	commonFaces.coordinates() = _coordinates;

	size_t count = 0;
	for (size_t i = 0; i < parts(); i++) {
		count += elementsCount[i];
	}
	commonFaces._elements.reserve(count);
	commonFaces._partPtrs.clear();
	commonFaces._partPtrs.reserve(_partPtrs.size());

	// create surface mesh
	commonFaces._partPtrs.push_back(commonFaces._elements.size());
	for (size_t i = 0; i + 1 < _partPtrs.size(); i++) {
		for (size_t j = 0; j < faces[i].size(); j++) {
			std::vector<eslocal> &face = faces[i][j];
			if (faces[i][j].size() == 3) {
				commonFaces._elements.push_back(new Triangle(&face[0]));
			}
			if (faces[i][j].size() == 4) {
				commonFaces._elements.push_back(new Square(&face[0]));
			}
		}
		commonFaces._partPtrs.push_back(commonFaces._elements.size());
		commonFaces.computeLocalIndices(commonFaces._partPtrs.size() - 1);
	}
}

void Mesh::getCornerLines(CornerLinesMesh &cornerLines) const
{
	if (parts() < 1) {
		std::cerr << "Internal error: _partPtrs.size()\n";
		exit(EXIT_FAILURE);
	}

	std::vector<std::vector<eslocal> > nodesElements(_coordinates.size());
	for (size_t i = 0; i < parts(); i++) {
		// Compute nodes' adjacent elements
		for (eslocal j = _partPtrs[i]; j < _partPtrs[i + 1]; j++) {
			for (size_t k = 0; k < _elements[j]->size(); k++) {
				nodesElements[_coordinates.clusterIndex(_elements[j]->node(k), i)].push_back(j);
			}
		}
	}

	std::vector<std::set<eslocal> > nodeParts(_coordinates.size());

	for (size_t i = 0; i < parts(); i++) {
		for (size_t j = 0; j < _coordinates.localSize(i); j++) {
			nodeParts[_coordinates.clusterIndex(j, i)].insert(i);
		}
	}

	for (size_t i = 0; i < parts(); i++) {
		for (eslocal j = _partPtrs[i]; j < _partPtrs[i + 1]; j++) {
			for (size_t k = 0; k < _elements[j]->faces(); k++) {
				std::vector<eslocal> face = _elements[j]->getFace(k);
				for (size_t f = 0; f < face.size(); f++) {
					face[f] = _coordinates.clusterIndex(face[f], i);
				}
				if (isOuterFace(nodesElements, face)) {
					for (size_t f = 0; f < face.size(); f++) {
						nodeParts[face[f]].insert(parts());
					}
				}
			}
		}
	}

	std::vector<std::set<eslocal> > neighbours(_coordinates.size());
	for (size_t p = 0; p < parts(); p++) {
		for (eslocal j = _partPtrs[p]; j < _partPtrs[p + 1]; j++) {
			for (size_t k = 0; k < _elements[j]->size(); k++) {
				std::vector<eslocal> n = _elements[j]->getNeighbours(k);
				eslocal c1 = _coordinates.clusterIndex(_elements[j]->node(k), p);
				for (size_t i = 0; i < n.size(); i++) {
					eslocal c2 = _coordinates.clusterIndex(n[i], p);
					if (c1 < c2) {
						neighbours[c1].insert(c2);
					}
				}
			}
		}
	}

	cornerLines.coordinates() = _coordinates;

	std::vector<eslocal> result;
	std::vector<eslocal>::iterator end;
	std::vector<eslocal>::iterator pit;
	std::vector<std::vector<eslocal> > pairs(parts());

	std::set<eslocal>::iterator it;
	for (size_t i = 0; i < neighbours.size(); i++) {
		for (it = neighbours[i].begin(); it != neighbours[i].end(); ++it) {
			result.resize(nodeParts[i].size());
			end = std::set_intersection(
				nodeParts[i].begin(), nodeParts[i].end(),
				nodeParts[*it].begin(), nodeParts[*it].end(),
				result.begin());
			if (end - result.begin() >= 3) {
				for (pit = result.begin(); pit != end && *pit < parts(); ++pit) {
					pairs[*pit].push_back(i);
					pairs[*pit].push_back(*it);
				}
			}
		}
	}

	eslocal size = 0;
	for (size_t i = 0; i < parts(); i++) {
		size += pairs[i].size() / 2;
	}

	cornerLines._elements.reserve(size);
	cornerLines._partPtrs.clear();
	cornerLines._partPtrs.reserve(_partPtrs.size());

	cornerLines._partPtrs.push_back(cornerLines._elements.size());
	for (size_t i = 0; i < parts(); i++) {
		for (size_t j = 0; j < pairs[i].size(); j += 2) {
			cornerLines._elements.push_back(new Line(&pairs[i][j]));
		}
		cornerLines._partPtrs.push_back(cornerLines._elements.size());
		cornerLines.computeLocalIndices(cornerLines._partPtrs.size() - 1);
	}

}

void Mesh::computeCorners(eslocal number, bool vertex, bool edges, bool faces)
{
	if (parts() < 1) {
		std::cerr << "Internal error: _partPtrs.size()\n";
		exit(EXIT_FAILURE);
	}
	if (parts() == 1 || (!vertex && !edges && !faces)) {
		return;
	}

	// node to element
	std::vector<std::vector<eslocal> > nodesElements(_coordinates.size());
	// node to neighbors nodes
	std::vector<std::set<eslocal> > neighbours(_coordinates.size());

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
	std::vector<eslocal> projection(_coordinates.size(), 0);

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
	Mesh cfm(_rank, _size);
	Element *el;
	eslocal index = 0;
	for (size_t i = 0; i < _coordinates.size(); i++) {
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
	//cfm.coordinates() = _coordinates;
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
	cfm._coordinates.localClear();
	cfm._coordinates.localResize(cfm.parts());
	for (size_t p = 0; p < cfm.parts(); p++) {
		cfm.computeLocalIndices(p);
	}

	if (faces) {
		cfm.computeFixPoints(number);
		for (size_t p = 0; p < cfm.parts(); p++) {
			for (size_t i = 0; i < cfm.getFixPoints()[p].size(); i++) {
				_subdomainBoundaries.setCorner(cfm.coordinates().globalIndex(cfm.getFixPoints()[p][i], p));
			}
		}
	}

	if (!edges && !vertex) {
		return;
	}

	Mesh clm(_rank, _size);
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

	std::vector<std::pair<eslocal, eslocal> > outerFaces;
	std::vector<char> nSubdomains(cfm.coordinates().size() * cfm.parts(), 0);
	std::vector<eslocal> points(cfm.coordinates().size(), 0);
	clm._partPtrs.clear();
	clm._partPtrs.push_back(0);
	for (size_t p = 0; p < cfm.parts(); p++) {
		for (size_t e = cfm._partPtrs[p]; e < cfm._partPtrs[p + 1]; e++) {
			for (size_t f = 0; f < cfm._elements[e]->faces(); f++) {
				std::vector<eslocal> face = cfm._elements[e]->getFace(f);
				if (isOuterFace(nodesFaces[p], face)) {
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

	eslocal linePoints = 0;
	for (size_t i = 0; i < cfm._coordinates.size(); i++) {
		if (points[i] == 1) {
			clm.coordinates().add(cfm._coordinates[i], linePoints, i);
			points[i] = linePoints++;
		}
	}

	std::vector<bool> mask(cfm.coordinates().size(), false);
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
	while (maskCounter < cfm.coordinates().size()) {
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
		for (size_t i = 0; i < outerFaces.size(); i++) {
			char *p0s = nSubdomains.data() + outerFaces[i].first * cfm.parts();
			char *p1s = nSubdomains.data() + outerFaces[i].second * cfm.parts();
			char *p0e = nSubdomains.data() + (outerFaces[i].first + 1) * cfm.parts();
			char *p1e = nSubdomains.data() + (outerFaces[i].second + 1) * cfm.parts();
			if (std::equal(s, e, p0s) && std::equal(s, e, p1s)) {
				tmpPair[0] = clm.coordinates().clusterIndex(outerFaces[i].first);
				tmpPair[1] = clm.coordinates().clusterIndex(outerFaces[i].second);
				clm._elements.push_back(new Line(tmpPair.data()));
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
		if (vertex) {
			for (auto it = vertices.begin(); it != vertices.end(); ++it) {
				eslocal fPoint = clm.coordinates().globalIndex(*it);
				_subdomainBoundaries.setCorner(cfm.coordinates().globalIndex(fPoint));
			}
		}
		if (clm._partPtrs.back() < clm._elements.size()) {
			clm._partPtrs.push_back(clm._elements.size());
		}
	}

	clm._coordinates.localClear();
	clm._coordinates.localResize(clm.parts());
	for (size_t p = 0; p < clm.parts(); p++) {
		clm.computeLocalIndices(p);
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

std::ostream& mesh::operator<<(std::ostream& os, const Mesh &m)
{
	for (size_t i = 0; i < m._elements.size(); i++) {
		os << *(m._elements[i]) << "\n";
	}
	return os;
}


