#include "mesh.h"

Mesh::Mesh(Coordinates &coordinates)
	:_coordinates(coordinates), _elements(0), _lastNode(0), _partsNodesCount(1, 0),
	 _fixPoints(0), _flags(flags::FLAGS_SIZE, false), _maxElementSize(0)
{
	_partPtrs.resize(2);
	_partPtrs[0] = 0;
	_partPtrs[1] = 0;
}

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

Mesh::Mesh(const Mesh &other)
	:_coordinates(other._coordinates), _lastNode(other._lastNode), _partPtrs(other._partPtrs),
	 _partsNodesCount(other._partsNodesCount), _fixPoints(other._fixPoints),
	 _flags(other._flags), _maxElementSize(other._maxElementSize)
{
	_elements.reserve(other._elements.size());
	for (size_t i = 0; i < other._elements.size(); i++) {
		_elements.push_back(other._elements[i]->copy());
	}
}

Mesh& Mesh::operator=(const Mesh &other)
{
	if (this != &other) {
		Mesh mesh(other);
		Mesh::assign(*this, mesh);
	}
	return *this;
}

void Mesh::assign(Mesh &m1, Mesh &m2)
{
	m1._elements.swap(m2._elements);
	m1._fixPoints.swap(m2._fixPoints);
	m1._flags.swap(m2._flags);
	m1._partPtrs.swap(m2._partPtrs);
	m1._partsNodesCount.swap(m2._partsNodesCount);
	m1._maxElementSize = m2._maxElementSize;
	m1._lastNode = m2._lastNode;
}

void Mesh::reserve(size_t size)
{
	_elements.reserve(size);
}

void Mesh::pushElement(Element* e)
{
	for (size_t i = 0; i < e->size(); i++) {
		if (_lastNode < e->node(i)) {
			_lastNode = e->node(i);
		}
	}
	if (_maxElementSize < e->size()) {
		_maxElementSize = e->size();
	}
	_elements.push_back(e);
	if (_flags[flags::NEW_PARTITION]) {
		_partPtrs.push_back(_partPtrs.back());
		_partsNodesCount.push_back(0);
		_flags[flags::NEW_PARTITION] = false;
	}
	_partPtrs.back() = _elements.size();
}

void Mesh::endPartition()
{
	computeLocalIndices(_partsNodesCount.size() - 1);
	_flags[flags::NEW_PARTITION] = true;
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
		std::cerr << "Unknown element with indices: ";
		for (idx_t i = 0; i < n; i++) {
			std::cerr << indices[i] << " ";
		}
		std::cerr << "\n";
		exit(EXIT_FAILURE);
	}

	if (_maxElementSize < e->size()) {
		_maxElementSize = e->size();
	}
	return e;
}

void Mesh::_assemble_matrix(SparseVVPMatrix &K, SparseVVPMatrix &M, std::vector<double> &f, idx_t part, bool dynamic)
{
	int nK = _partsNodesCount[part] * Point::size();
	K.resize(nK, nK);
	f.resize(nK);

	size_t maxElementSize = _maxElementSize * Point::size();
	std::vector<double> Ke(maxElementSize * maxElementSize);
	std::vector<double> fe(maxElementSize);

	std::vector <double> inertia (3, 0.0);
	inertia[2] = 9810.0 * 7.85e-9;
	double ex = 2.1e5;
	double mi = 0.3;

	if (dynamic) {
		M.resize(nK, nK);
		// Only one block is assembled -> all others are the same
		std::vector <double> Me(_maxElementSize * _maxElementSize);

		for (int i = _partPtrs[part]; i < _partPtrs[part + 1]; i++) {
			_elements[i]->elasticity(Ke, Me, fe, _coordinates, inertia, ex, mi);
			_elements[i]->addLocalValues(K, M, f, Ke, Me, fe, 0);
		}
	} else {
		for (int i = _partPtrs[part]; i < _partPtrs[part + 1]; i++) {
			_elements[i]->elasticity(Ke, fe, _coordinates, inertia, ex, mi);
			_elements[i]->addLocalValues(K, f, Ke, fe, 0);
		}
	}
}

void Mesh::partitiate(idx_t parts, idx_t fixPoints)
{
	_partPtrs.resize(parts + 1);
	_partsNodesCount.resize(parts);

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

	cilk_for (idx_t i = 0; i < parts; i++) {
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
		computeLocalIndices(part);
	}
}

void Mesh::computeLocalIndices(size_t part)
{
	idx_t *nodeMap = new idx_t[_lastNode + 1];

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


void Mesh::saveBasis(std::ofstream &vtk, std::vector<std::vector<int> > &l2g_vec)
{
	vtk.open("mesh.vtk", std::ios::out | std::ios::trunc);
	vtk << "# vtk DataFile Version 3.0\n";
	vtk << "Test\n";
	vtk << "ASCII\n\n";
	vtk << "DATASET UNSTRUCTURED_GRID\n";
	size_t nSubClst = l2g_vec.size();
	size_t cnt=0;


	size_t n_points = 0;
	for (size_t d = 0; d < l2g_vec.size(); d++) {
		n_points+=l2g_vec[d].size();
	}

	vtk << "POINTS " << n_points << " float\n";
	for (size_t d = 0; d < nSubClst; d++) {
		for (size_t i = 0; i < l2g_vec[d].size(); i++) {
			vtk << _coordinates[l2g_vec[d][i]].x << " " ;
			vtk << _coordinates[l2g_vec[d][i]].y << " " ;
			vtk << _coordinates[l2g_vec[d][i]].z << "\n";
		}
	}

	size_t size = 0;
	for (size_t i = 0; i < _elements.size(); i++) {
		size += _elements[i]->size() + 1;
	}
	vtk << "CELLS " << _elements.size() << " " << size << "\n";

	size_t i=0;
	for (size_t part = 0; part + 1 < _partPtrs.size(); part++) {
		for (idx_t ii = 0; ii < _partPtrs[part + 1] - _partPtrs[part]; ii++) {
			vtk << _elements[i]->size();
			for (size_t j = 0; j < _elements[i]->size(); j++) {
				vtk << " " << _elements[i]->localNode(j)+ cnt;// - _coordinates.getOffset() ;
			}
			vtk << "\n";
			i++;
		}
		cnt += l2g_vec[part].size();
	}

	vtk << "\n";
	vtk << "CELL_TYPES " << _elements.size() << "\n";
	for (size_t i = 0; i < _elements.size(); i++) {
		vtk << _elements[i]->vtkCode() << "\n";
	}

	vtk << "\n";
	vtk << "CELL_DATA " << _elements.size() << "\n";
	vtk << "SCALARS decomposition int 1\n";
	vtk << "LOOKUP_TABLE decomposition\n";
	for (size_t part = 0; part + 1 < _partPtrs.size(); part++) {
		for (idx_t i = 0; i < _partPtrs[part + 1] - _partPtrs[part]; i++) {
			vtk << part << "\n";
		}
	}
}

void Mesh::saveVTK(std::vector<std::vector<double> > &displacement, std::vector<std::vector <int> > &l2g_vec)
{
	std::ofstream vtk;
	saveBasis(vtk, l2g_vec);

	size_t n_points = 0;
	for (size_t d = 0; d < l2g_vec.size(); d++) {
		n_points += l2g_vec[d].size();
	}

	vtk << "\n";
	vtk << "POINT_DATA " << n_points << "\n";
	vtk << "SCALARS displacements float 3\n";
	vtk << "LOOKUP_TABLE default\n";
	for (size_t i = 0; i < displacement.size(); i++) {
		for (size_t j = 0; j < displacement[i].size() / 3; j++) {
			vtk << displacement[i][3 * j + 0] << " ";
			vtk << displacement[i][3 * j + 1] << " ";
			vtk << displacement[i][3 * j + 2] << "\n";
		}
		
	}

	vtk.close();
}

void Mesh::saveVTK(std::vector<std::vector <int> > &l2g_vec)
{
	std::ofstream vtk;
	saveBasis(vtk, l2g_vec);

	vtk.close();
}

std::ostream& operator<<(std::ostream& os, const Mesh &m)
{
	for (size_t i = 0; i < m._elements.size(); i++) {
		os << *(m._elements[i]) << "\n";
	}
	return os;
}
