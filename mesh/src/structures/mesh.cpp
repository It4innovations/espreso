#include "mesh.h"

Mesh::Mesh(Coordinates &coordinates)
	:_coordinates(coordinates), _indicesType(Element::GLOBAL), _elements(0), _lastNode(0),
	 _partsNodesCount(1, 0), _fixPoints(0), _flags(flags::FLAGS_SIZE, false), _maxElementSize(0)
{
	_partPtrs.resize(2);
	_partPtrs[0] = 0;
	_partPtrs[1] = 0;
}

Mesh::Mesh(const char *fileName, Coordinates &coordinates, idx_t parts, idx_t fixPoints):
	_indicesType(Element::GLOBAL),
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
	:_indicesType(other._indicesType), _coordinates(other._coordinates),
	 _lastNode(other._lastNode), _partPtrs(other._partPtrs),
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
	m1._coordinates = m2._coordinates;
	m1._indicesType = m2._indicesType;
	m1._elements.swap(m2._elements);
	m1._lastNode = m2._lastNode;
	m1._partPtrs.swap(m2._partPtrs);
	m1._partsNodesCount.swap(m2._partsNodesCount);
	m1._fixPoints.swap(m2._fixPoints);
	m1._flags.swap(m2._flags);
	m1._maxElementSize = m2._maxElementSize;
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
	if (Tetrahedron4::match(indices, n)) {
		e = new Tetrahedron4(indices);
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

void Mesh::_elasticity(SparseVVPMatrix &K, SparseVVPMatrix &M, std::vector<double> &f, idx_t part, bool dynamic)
{
	int nK = _partsNodesCount[part] * Point::size();
	K.resize(nK, nK);
	f.resize(nK);

	size_t maxElementSize = _maxElementSize * Point::size();
	std::vector<double> Ke(maxElementSize * maxElementSize);
	std::vector<double> Me;
	std::vector<double> fe(maxElementSize);

	std::vector <double> inertia (3, 0.0);
	inertia[2] = 9810.0 * 7.85e-9;
	double ex = 2.1e5;
	double mi = 0.3;

	if (dynamic) {
		M.resize(nK, nK);
		// Only one block is assembled -> all others are the same
		Me.resize(_maxElementSize * _maxElementSize);
	}

	for (int i = _partPtrs[part]; i < _partPtrs[part + 1]; i++) {
		_assembleElesticity(_elements[i], part, Ke, Me, fe, inertia, ex, mi, dynamic);
		_integrateElasticity(_elements[i], K, M, f, Ke, Me, fe, dynamic);
	}
}

void Mesh::_assembleElesticity(
		const Element *e,
		size_t part,
		std::vector<double> &Ke,
		std::vector<double> &Me,
		std::vector<double> &fe,
		std::vector<double> &inertia,
		double ex,
		double mi,
		bool dynamic) const
{
	const std::vector<std::vector<double> > &dN = e->dN();
	const std::vector<std::vector<double> > &N = e->N();
	const std::vector<double> &weighFactor = e->weighFactor();

	std::vector<double> coordinates(Point::size() * e->size());
	const std::vector<idx_t> &l2g = _coordinates.localToGlobal(part);
	for (size_t i = 0; i < e->size(); i++) {
		&coordinates[i * Point::size()] << _coordinates[l2g[e->node(i)]];
	}

	int dimension = Point::size();
	int Ksize = dimension * e->size();
	int Csize = 6;	// TODO: even for D2??
	int nodes = e->size();
	int gausePoints = e->gpSize();

	std::vector<double> MatC(Csize * Csize, 0.0);

	double E = ex / ((1 + mi) * (1 - 2 * mi));
	double mi2 = E * (1.0 - mi);
	double mi3 = E * (0.5 - mi);
	mi = E * mi;

	MatC[0]  = mi2; MatC[1]  = mi;  MatC[2]  = mi;
	MatC[6]  = mi;  MatC[7]  = mi2; MatC[8]  = mi;
	MatC[12] = mi;  MatC[13] = mi;  MatC[14] = mi2;

	MatC[21] = mi3; MatC[28] = mi3; MatC[35] = mi3;


	std::vector<double> MatJ (dimension * dimension, 0);
	std::vector<double> invJ (dimension * dimension, 0);

	std::vector<double> dND (Ksize, 0);
	std::vector<double> CB (Ksize * Csize, 0);
	std::vector<double> B (Ksize * Csize, 0);

	Ke.resize(Ksize * Ksize);
	fe.resize(Ksize);
	fill(Ke.begin(), Ke.end(), 0);
	fill(fe.begin(), fe.end(), 0);
	if (dynamic) {
		Me.resize(nodes * nodes);
		fill(Me.begin(), Me.end(), 0);
	}

	for (int gp = 0; gp < gausePoints; gp++) {

		cblas_dgemm(
			CblasRowMajor, CblasNoTrans, CblasNoTrans,
			dimension, dimension, nodes,
			1, &dN[gp][0], nodes, &coordinates[0], dimension,
			0, &MatJ[0], dimension
		);

		double detJ = fabs( MatJ[0] * MatJ[4] * MatJ[8] +
							MatJ[1] * MatJ[5] * MatJ[6] +
							MatJ[2] * MatJ[3] * MatJ[7] -
							MatJ[2] * MatJ[4] * MatJ[6] -
							MatJ[1] * MatJ[3] * MatJ[8] -
							MatJ[0] * MatJ[5] * MatJ[7]);

		double detJx = 1 / detJ;

		invJ[0] = detJx * (  MatJ[8] * MatJ[4] - MatJ[7] * MatJ[5] );
		invJ[1] = detJx * (- MatJ[8] * MatJ[1] + MatJ[7] * MatJ[2] );
		invJ[2] = detJx * (  MatJ[5] * MatJ[1] - MatJ[4] * MatJ[2] );
		invJ[3] = detJx * (- MatJ[8] * MatJ[3] + MatJ[6] * MatJ[5] );
		invJ[4] = detJx * (  MatJ[8] * MatJ[0] - MatJ[6] * MatJ[2] );
		invJ[5] = detJx * (- MatJ[5] * MatJ[0] + MatJ[3] * MatJ[2] );
		invJ[6] = detJx * (  MatJ[7] * MatJ[3] - MatJ[6] * MatJ[4] );
		invJ[7] = detJx * (- MatJ[7] * MatJ[0] + MatJ[6] * MatJ[1] );
		invJ[8] = detJx * (  MatJ[4] * MatJ[0] - MatJ[3] * MatJ[1] );

		cblas_dgemm(
			CblasRowMajor, CblasNoTrans, CblasNoTrans,
			dimension, nodes, dimension,
			1, &invJ[0], dimension, &dN[gp][0], nodes,
			0, &dND.front(), nodes
		);

		// TODO: block ordering inside B
		int columns = Ksize;
		const double *dNDx = &dND[0];
		const double *dNDy = &dND[e->size()];
		const double *dNDz = &dND[2 * e->size()];

		// B =
		// dX   0   0
		//  0  dY   0
		//  0   0  dZ
		// dY  dX   0
		//  0  dZ  dY
		// dZ   0  dX

		memcpy(&B[0],                           dNDx, sizeof(double) * e->size());
		memcpy(&B[3 * columns + e->size()],     dNDx, sizeof(double) * e->size());
		memcpy(&B[5 * columns + 2 * e->size()], dNDx, sizeof(double) * e->size());

		memcpy(&B[1 * columns + e->size()],     dNDy, sizeof(double) * e->size());
		memcpy(&B[3 * columns],                 dNDy, sizeof(double) * e->size());
		memcpy(&B[4 * columns + 2 * e->size()], dNDy, sizeof(double) * e->size());

		memcpy(&B[2 * columns + 2 * e->size()], dNDz, sizeof(double) * e->size());
		memcpy(&B[4 * columns + e->size()],     dNDz, sizeof(double) * e->size());
		memcpy(&B[5 * columns],                 dNDz, sizeof(double) * e->size());


		// C * B
		cblas_dgemm(
			CblasRowMajor, CblasNoTrans, CblasNoTrans,
			Csize, Ksize, Csize,
			1, &MatC[0], Csize, &B[0], Ksize,
			0, &CB.front(), Ksize
		);
		//Ke = Ke + (B' * (C * B)) * dJ * WF;
		cblas_dgemm(
			CblasRowMajor, CblasTrans, CblasNoTrans,
			Ksize, Ksize, Csize,
			detJ * weighFactor[gp], &B[0], Ksize, &CB[0], Ksize,
			1, &Ke[0], Ksize
		);

		for (int i = 0; i < Ksize; i++) {
			fe[i] += detJ * weighFactor[gp] * N[gp][i % nodes] * inertia[i / nodes];
		}

		if (dynamic) {
			// Me = Me + WF * (DENS * dJ) * (N' * N);
			double dense = 7.85e-9;
			cblas_dgemm(
				CblasRowMajor, CblasTrans, CblasNoTrans,
				nodes, nodes, 1,
				dense * detJ * weighFactor[gp], &N[gp][0], nodes, &N[gp][0], nodes,
				1, &Me[0], nodes
			);
		}
	}
}

void Mesh::_integrateElasticity(
		const Element *e,
		SparseVVPMatrix &K,
		SparseVVPMatrix &M,
		std::vector<double> &f,
		const std::vector<double> &Ke,
		const std::vector<double> &Me,
		const std::vector<double> &fe,
		bool dynamic
	) const
{
	// Element ordering: xxxx, yyyy, zzzz,...
	// Global ordering:  xyz, xyz, xyz, xyz, ...
	size_t row, column;
	size_t s = Point::size();

	for (size_t i = 0; i < s * e->size(); i++) {
		row = s * (e->node(i % e->size())) + i / e->size();
		for (size_t j = 0; j < s * e->size(); j++) {
			column = s * (e->node(j % e->size())) + j / e->size();
			K(row, column) = Ke[i * s * e->size() + j];
		}
		f[row] += fe[i];
	}
	if (!dynamic) {
		return;
	}
	for (size_t i = 0; i < e->size(); i++) {
		row = s * (e->node(i));
		for (size_t j = 0; j < e->size(); j++) {
			column = s * (e->node(i));
			for (size_t k = 0; k < s; k++) {
				M(row + k, column + k) += Me[i * e->size() + j];
			}
		}
	}
}

void Mesh::partitiate(idx_t parts, idx_t fixPoints)
{
	_partPtrs.resize(parts + 1);
	_partsNodesCount.resize(parts);
	_coordinates.localClear();
	_coordinates.localResize(parts);

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

#ifndef SEQUENTIAL
	cilk_for (idx_t i = 0; i < parts; i++) {
#else
	for (idx_t i = 0; i < parts; i++) {
#endif
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
		_elements[i]->fillNodes(n + e[index]);
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
	_indicesType = Element::LOCAL;
}

void Mesh::computeLocalIndices(size_t part)
{
	std::vector<idx_t> nodeMap (_lastNode + 1, -1);

	// Compute mask of nodes
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

	_coordinates.computeLocal(part, nodeMap, nSize);
}

idx_t Mesh::getCentralNode(idx_t first, idx_t last, idx_t *ePartition, idx_t part, idx_t subpart) const
{
	// Compute CSR format of symmetric adjacency matrix
	////////////////////////////////////////////////////////////////////////////
	std::vector<std::set<int> > neighbours(_partsNodesCount[part]);
	for (idx_t i = first, index = 0; i < last; i++, index++) {
		if (ePartition[index] == subpart) {
			for (size_t j = 0; j < _elements[i]->size(); j++) {
				std::vector<idx_t> neigh = _elements[i]->getNeighbours(j);
				for (size_t k = 0; k < neigh.size(); k++) {
					if (_elements[i]->node(j) < neigh[k]) {
						neighbours[_elements[i]->node(j)].insert(neigh[k]);
					}
				}
			}
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

void Mesh::saveNodeArray(idx_t *nodeArray, size_t part)
{
	size_t p = 0;
	for (idx_t i = _partPtrs[part]; i < _partPtrs[part + 1]; i++) {
		_elements[i]->fillNodes(nodeArray + p);
		p += _elements[i]->size();
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
	size_t cnt = 0;


	size_t n_points = 0;
	for (size_t d = 0; d < l2g_vec.size(); d++) {
		n_points += l2g_vec[d].size();
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
				vtk << " " << _elements[i]->node(j)+ cnt;
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

void Mesh::saveVTK(const char* filename)
{
	std::ofstream vtk;

	vtk.open(filename, std::ios::out | std::ios::trunc);
	vtk << "# vtk DataFile Version 3.0\n";
	vtk << "Test\n";
	vtk << "ASCII\n\n";
	vtk << "DATASET UNSTRUCTURED_GRID\n";
	vtk << "POINTS " << _coordinates.size() << " float\n";
	vtk << _coordinates << "\n";

	size_t size = 0;
	for (size_t i = 0; i < _elements.size(); i++) {
		size += _elements[i]->size() + 1;
	}

	vtk << "CELLS " << _elements.size() << " " << size << "\n";
	for (size_t p = 0; p + 1 < _partPtrs.size(); p++) {
		std::vector<idx_t> l2g = _coordinates.localToGlobal(p);
		for (size_t i = _partPtrs[p]; i < _partPtrs[p + 1]; i++) {
			vtk << _elements[i]->size();
			for (size_t j = 0; j < _elements[i]->size(); j++) {
				vtk << " " << l2g[_elements[i]->node(j)] - _coordinates.getOffset();
			}
			vtk << "\n";
		}
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

	vtk.close();
}

bool Mesh::isOuterFace(std::vector<std::vector<int> > &nodesElements, std::vector<idx_t> &face)
{
	std::vector<int> result(nodesElements[face[0]]);
	std::vector<int>::iterator it = result.end();

	for (size_t i = 1; i < face.size(); i++) {
		std::vector<int> tmp(result.begin(), it);
		it = std::set_intersection(
		         tmp.begin(), tmp.end(),
		         nodesElements[face[i]].begin(), nodesElements[face[i]].end(),
		         result.begin());
		if (it - result.begin() == 1) {
			return true;
		}
	}
	return false;
}

void Mesh::getBoundary(BoundaryMesh &boundaryMesh)
{
	std::vector<std::vector<std::vector<idx_t> > > faces(_partPtrs.size() - 1);
	std::vector<size_t> elementsCount(_partPtrs.size() - 1, 0);
	std::vector<idx_t> selection(_coordinates.size() + _coordinates.getOffset(), -1);

	if (_partPtrs.size() < 2) {
		std::cerr << "Internal error: _partPtrs.size()\n";
		exit(EXIT_FAILURE);
	}
#ifndef SEQUENTIAL
	cilk_for (size_t i = 0; i < _partPtrs.size() - 1; i++) {
#else
	for (size_t i = 0; i < _partPtrs.size() - 1; i++) {
#endif
		// Compute nodes' adjacent elements
		const std::vector<idx_t> &l2g = _coordinates.localToGlobal(i);
		std::vector<std::vector<int> > nodesElements(l2g.size());
		for (idx_t j = _partPtrs[i]; j < _partPtrs[i + 1]; j++) {
			for (size_t k = 0; k < _elements[j]->size(); k++) {
				nodesElements[_elements[j]->node(k)].push_back(j);
			}
		}


		for (idx_t j = _partPtrs[i]; j < _partPtrs[i + 1]; j++) {
			for (size_t k = 0; k < _elements[j]->faces(); k++) {
				std::vector<idx_t> face = _elements[j]->getFace(k);
				if (isOuterFace(nodesElements, face)) {
					faces[i].push_back(face);
					if (face.size() == 3) {
						elementsCount[i] += 1;
					}
					if (face.size() == 4) {
						elementsCount[i] += 2;
					}
					for (size_t p = 0; p < face.size(); p++) {
						selection[l2g[face[p]]] = 1;
					}
				}
			}
		}
	}

	Coordinates &coords = boundaryMesh.coordinates();
	coords.localResize(_partPtrs.size() - 1);

	size_t c = 0;
	coords.resize(std::count(selection.begin(), selection.end(), 1));
	for (size_t i = 0; i < selection.size(); i++) {
		if (selection[i] == 1) {
			selection[i] = c;
			coords[c++] = _coordinates[i];
		}
	}

#ifndef SEQUENTIAL
	cilk_for (size_t i = 0; i < _partPtrs.size() - 1; i++) {
#else
	for (size_t i = 0; i < _partPtrs.size() - 1; i++) {
#endif
		const std::vector<idx_t> &l2g = _coordinates.localToGlobal(i);
		for (size_t j = 0; j < faces[i].size(); j++) {
			for (size_t k = 0; k < faces[i][j].size(); k++) {
				faces[i][j][k] = selection[l2g[faces[i][j][k]]];
			}
		}
	}

	size_t count = 0;
	for (size_t i = 0; i + 1 < _partPtrs.size(); i++) {
		count += elementsCount[i];
	}
	boundaryMesh.reserve(count);

	for (size_t i = 0; i + 1 < _partPtrs.size(); i++) {
		for (size_t j = 0; j < faces[i].size(); j++) {
			std::vector<idx_t> &face = faces[i][j];
			if (faces[i][j].size() == 3) {
				boundaryMesh.pushElement(new Triangle(&face[0]));
			}
			if (faces[i][j].size() == 4) {
				size_t min = 0;
				for (size_t p = 1; p < 4; p++) {
					if (face[min] > face[p]) {
						min = p;
					}
				}
				if (min % 2 == 0) {
					boundaryMesh.pushElement(new Triangle(&face[0]));
					face[1] = face[0];
					boundaryMesh.pushElement(new Triangle(&face[1]));
				} else {
					boundaryMesh.pushElement(new Triangle(&face[1]));
					face[2] = face[3];
					boundaryMesh.pushElement(new Triangle(&face[0]));
				}
			}
		}
		boundaryMesh.endPartition();
	}
}

void BoundaryMesh::elasticity(DenseMatrix &K, size_t part) const
{
	idx_t nK = Point::size() * _partsNodesCount[part];
	int eSize = _partPtrs[part + 1] - _partPtrs[part];
	K.resize(nK, nK);
	std::vector<double> nodes(nK);
	std::vector<int> elems(3 * eSize);

	const std::vector<idx_t> &l2g = _coordinates.localToGlobal(part);
	for (size_t i = 0; i < _partsNodesCount[part]; i++) {
		&nodes[i * Point::size()] << _coordinates[l2g[i]];
	}
	for (size_t i = _partPtrs[part]; i < _partPtrs[part + 1]; i++) {
		// TODO: various data types int32_t and int64_t
		// _elements[i]->fillNodes(&elems[3 * i]); CANNOT be used
		for (size_t j = 0; j < _elements[i]->size(); j++) {
			elems[3 * i + j] = _elements[i]->node(j);
		}
	}

	/*bem4i::getLameSteklovPoincare(
	    K.values(),
	    _partsNodesCount[part],
	    &nodes[0],
	    eSize,
	    &elems[0],
	    0.29,			// nu
	    1e11,			// E
	    3,				// order near
	    4,				// order far
	    true			// verbose
	    );*/
}

std::ostream& operator<<(std::ostream& os, const Mesh &m)
{
	for (size_t i = 0; i < m._elements.size(); i++) {
		os << *(m._elements[i]) << "\n";
	}
	return os;
}
