#include "mesh.h"

Mesh::Mesh(const char *fileName, Coordinates &coordinates, idx_t parts, idx_t fixPoints):
	_coordinates(coordinates),
	_flags(flags::FLAGS_SIZE, false),
	_maxElementSize(0)
{
	_elements.resize(Loader::getLinesCount(fileName));

	std::ifstream file(fileName);
	std::string line;

	idx_t indices[20], n; 	// 20 is the max of vertices of a element
	double value;
	idx_t minIndices = 10000;

	if (file.is_open()) {
		_lastNode = 0;
		for (idx_t c = 0; c < _elements.size(); c++) {
			getline(file, line, '\n');;
			std::stringstream ss(line);
			n = 0;
			while(ss >> value) {
				indices[n++] = value;
				if (minIndices > value) {
					minIndices = value;
				}
				if (_lastNode < value) {
					_lastNode = value;
				}
			}
			_elements[c] = createElement(indices, n);
		}
		file.close();
	} else {
		fprintf(stderr, "Cannot load mesh from file: %s.\n", fileName);
		exit(EXIT_FAILURE);
	}

	// correct indexing -> C/C++ indexes start at 0, but Points usually start at 1
	coordinates.setOffset(minIndices);

	if (parts > 0) {
		partitiate(parts, fixPoints);
	} else if (fixPoints > 0) {
		computeFixPoints(fixPoints);
	}
}

Element* Mesh::createElement(idx_t *indices, idx_t n)
{
	Element *e = NULL;
	if (Tetrahedron::match(indices, n)) {
		e = new Tetrahedron(indices);
	}
	if (Hexahedron::match(indices, n)) {
		e = new Hexahedron(indices);
	}

	if (e == NULL) {
		fprintf(stderr, "Unknown element with indices: ");
		for (idx_t i = 0; i < n; i++) {
			fprintf(stderr, "%ld ", indices[i]);
		}
		fprintf(stderr, ".\n");
		exit(EXIT_FAILURE);
	}

	if (_maxElementSize < e->size()) {
		_maxElementSize = e->size();
	}
	return e;
}

void Mesh::_assemble_matrix(SparseDOKMatrix &K, SparseDOKMatrix &M, std::vector<double> &f, idx_t part, bool dynamic)
{
	int nK = _partsNodesCount[part] * Point::size();
	SparseDOKMatrix DOK_K(nK, nK);	// temporary matrix used for construct K
	f.resize(nK);

	size_t maxElementSize = _maxElementSize * Point::size();
	std::vector<double> Ke(maxElementSize * maxElementSize);
	std::vector<double> fe(maxElementSize);

	std::vector <double> inertia (3, 0.0);
	inertia[2] = 9810.0 * 7.85e-9;
	double ex = 2.1e5;
	double mi = 0.3;

	if (dynamic) {
		SparseDOKMatrix DOK_M(nK, nK);	// temporary matrix used for construct M
		// Only one block is assembled -> all others are the same
		std::vector <double> Me(_maxElementSize * _maxElementSize);

		for (int i = _partPtrs[part]; i < _partPtrs[part + 1]; i++) {
			_elements[i]->elasticity(Ke, Me, fe, _coordinates, inertia, ex, mi);
			_elements[i]->addLocalValues(DOK_K, DOK_M, f, Ke, Me, fe, 0);
		}
		M = DOK_M;
	} else {
		for (int i = _partPtrs[part]; i < _partPtrs[part + 1]; i++) {
			_elements[i]->elasticity(Ke, fe, _coordinates, inertia, ex, mi);
			_elements[i]->addLocalValues(DOK_K, f, Ke, fe, 0);
		}
	}

	K = DOK_K;
}

void Mesh::partitiate(idx_t parts, idx_t fixPoints)
{
	_partPtrs.resize(parts + 1);
	_partsNodesCount.resize(parts + 1);

	// Call METIS to get partition of a whole mesh
	idx_t *elementPartition = getPartition(0, _elements.size(), parts);

	// Rearrange mesh's elements
	partitiate(elementPartition);

	delete[] elementPartition;

	if (fixPoints > 0) {
		computeFixPoints(fixPoints);
	}

	_flags[flags::PARTITIONS] = true;
}

void Mesh::computeFixPoints(idx_t fixPoints)
{
	idx_t parts = _partPtrs.size() - 1;
	_fixPoints.resize(parts * fixPoints);

	for (idx_t i = 0; i < parts; i++) {
		idx_t *eSubPartition = getPartition(_partPtrs[i], _partPtrs[i + 1], fixPoints);

		for (idx_t j = 0; j < fixPoints; j++) {
			_fixPoints[i * fixPoints + j] = getCentralNode(_partPtrs[i], _partPtrs[i + 1], eSubPartition, i, j);
		}

		delete[] eSubPartition;
	}

	_flags[flags::FIX_POINTS] = true;
}

idx_t* Mesh::getPartition(idx_t first, idx_t last, idx_t parts) const
{
	// INPUTS
	idx_t ncommon, eSize, nSize, *e, *n, options[METIS_NOPTIONS];

	// OUTPUTS
	idx_t objval, *ePartition, *nPartition;


	// FILL INPUT VARIABLES
	////////////////////////////////////////////////////////////////////////////

	// number of common nodes to be neighbor
	ncommon = 2;

	// There is probably BUG in METIS numbering or I do not understand documentation.
	// The solution is increase the size of 'nodesCount' and keep the default numbering
	//options[METIS_OPTION_NUMBERING] = coordinates.getNumbering();
	METIS_SetDefaultOptions(options);

	eSize = last - first;
	nSize = _lastNode + 1;

	// create array storing pointers to elements' nodes
	e = new idx_t[eSize + 1];
	e[0] = 0;
	for (idx_t i = first, index = 0; i < last; i++, index++) {
		e[index + 1] = e[index] + _elements[i]->size();
	}

	// create array of nodes
	n = new idx_t[e[eSize]];
	for (idx_t i = first, index = 0; i < last; i++, index++) {
		_elements[i]->fillNodes(n + e[index], Element::GLOBAL);
	}

	// PREPARE OUTPUT VARIABLES
	////////////////////////////////////////////////////////////////////////////
	ePartition = new idx_t[eSize];
	nPartition = new idx_t[nSize];

	int result = METIS_PartMeshDual(
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
	                 nPartition
	             );
	checkMETISResult(result);

	delete[] e;
	delete[] n;
	delete[] nPartition;

	return ePartition;
}

void Mesh::partitiate(idx_t *ePartition)
{
	_partPtrs[0] = 0;

	Element *e;
	idx_t p;
	for (size_t part = 0; part < _partPtrs.size() - 1; part++) {
		idx_t index = _partPtrs[part];	// index of last ordered element
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

	computeLocalIndices();
}

void Mesh::computeLocalIndices()
{
	idx_t *nodeMap = new idx_t[_lastNode + 1];

	for (idx_t part = 0; part < _partPtrs.size() - 1; part++) {
		// Compute mask of nodes
		memset(nodeMap, 0, (_lastNode + 1) * sizeof(idx_t));
		for (idx_t e = _partPtrs[part]; e < _partPtrs[part + 1]; e++) {
			for (size_t n = 0; n < _elements[e]->size(); n++) {
				nodeMap[_elements[e]->node(n)] = 1;
			}
		}

		// re-index nodes
		idx_t nSize = 0;
		for (idx_t k = 0; k <= _lastNode; k++) {
			if (nodeMap[k] == 1) {
				nodeMap[k] = nSize++;
			}
		}
		_partsNodesCount[part] = nSize;

		for (idx_t e = _partPtrs[part]; e < _partPtrs[part + 1]; e++) {
			_elements[e]->setLocalIndices(nodeMap);
		}
	}

	delete[] nodeMap;
}

idx_t Mesh::getCentralNode(idx_t first, idx_t last, idx_t *ePartition, idx_t part, idx_t subpart) const
{
	// Compute CSR format of symmetric adjacency matrix
	////////////////////////////////////////////////////////////////////////////
	BoundaryNodes neighbours(_partsNodesCount[part]);
	for (idx_t i = first, index = 0; i < last; i++, index++) {
		if (ePartition[index] == subpart) {
			_elements[i]->fillNeighbour(neighbours, Element::LOCAL);
		}
	}
	MKL_INT nonZeroValues = 0;
	for (size_t i = 0; i < neighbours.size(); i++) {
		nonZeroValues += neighbours[i].size();
	}

	float *a;
	MKL_INT *ia, *ja;
	a = new float[nonZeroValues];
	ia = new MKL_INT[_partsNodesCount[part] + 1];
	ja = new MKL_INT[nonZeroValues];

	MKL_INT i = 0;
	std::set<int>::const_iterator it;
	for (size_t n = 0; n < neighbours.size(); n++) {
		std::set<int> &values = neighbours[n];
		ia[n] = i;
		for (it = values.begin(); it != values.end(); ++it, i++) {
			a[i] = 1;
			ja[i] = *it;
		}
	}
	ia[neighbours.size()] = i;

	MKL_INT nSize = _partsNodesCount[part];
	float *x, *y, *swap;
	x = new float[nSize];
	y = new float[nSize];

	// Initial vector
	for (MKL_INT xi = 0; xi < nSize; xi++) {
		x[xi] = 1. / nSize;
	}
	float last_l = nSize, l = 1;
	MKL_INT incr = 1;

	while(fabs((l - last_l) / l) > 1e-6) {
		mkl_cspblas_scsrsymv("U", &nSize, a, ia, ja, x, y);
		last_l = l;
		l = snrm2(&nSize, y, &incr);
		cblas_sscal(nSize, 1 / l, y, incr);
		swap = x;
		x = y;
		y = swap;
	}
	idx_t result = cblas_isamax(nSize, x, incr);

	delete[] a;
	delete[] ia;
	delete[] ja;
	delete[] x;
	delete[] y;

	return result;
}

void Mesh::checkMETISResult(int result) const
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

void Mesh::checkMKLResult(MKL_INT result) const
{
	switch (result) {
	case 0:
		return;
	default:
		fprintf(stderr, "MKL error: %i.\n", result);
		exit(EXIT_FAILURE);
	}
}

Mesh::~Mesh()
{
	for (size_t i = 0; i < _elements.size(); i++) {
		delete _elements[i];
	}
}

std::ostream& operator<<(std::ostream& os, const Mesh &m)
{
	for (size_t i = 0; i < m._elements.size(); i++) {
		os << *(m._elements[i]) << "\n";
	}
	return os;
}
