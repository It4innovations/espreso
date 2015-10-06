#include "mesh.h"

using namespace mesh;

Mesh::Mesh():_elements(0), _fixPoints(0)
{
	_partPtrs.resize(2);
	_partPtrs[0] = 0;
	_partPtrs[1] = 0;
}

Mesh::Mesh(const Mesh &other):
	_coordinates(other._coordinates), _partPtrs(other._partPtrs),
	_fixPoints(other._fixPoints)
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
	m1._elements.swap(m2._elements);
	m1._partPtrs.swap(m2._partPtrs);
	m1._fixPoints.swap(m2._fixPoints);
}


void Mesh::_elasticity(
		SparseVVPMatrix<eslocal> &K,
		SparseVVPMatrix<eslocal> &M,
		std::vector<double> &f,
		eslocal part,
		bool dynamic) const
{
	eslocal nK = _coordinates.localSize(part) * Point::size();
	K.resize(nK, nK);
	if (dynamic) {
		M.resize(nK, nK);
	}
	f.resize(nK);

	DenseMatrix Ke, Me;
	std::vector<double> fe;

	std::vector<double> inertia(3, 0.0);
	inertia[2] = 9810.0 * 7.85e-9;
	double ex = 2.1e5;
	double mi = 0.3;
	double E = ex / ((1 + mi) * (1 - 2 * mi));
	DenseMatrix C(6, 6);

	C(0, 1) = C(0, 2) = C(1, 0) = C(1, 2) = C(2, 0) = C(2, 1) = E * mi;
	C(0, 0) = C(1, 1) = C(2, 2) = E * (1.0 - mi);
	C(3, 3) = C(4, 4) = C(5, 5) = E * (0.5 - mi);

	for (eslocal i = _partPtrs[part]; i < _partPtrs[part + 1]; i++) {
		_assembleElesticity(_elements[i], part, Ke, Me, fe, inertia, C, dynamic);
		_integrateElasticity(_elements[i], K, M, f, Ke, Me, fe, dynamic);
	}
}

inline double determinant3x3(DenseMatrix &m)
{
	const double *values = m.values();
	return fabs(
			values[0] * values[4] * values[8]
					+ values[1] * values[5] * values[6]
					+ values[2] * values[3] * values[7]
					- values[2] * values[4] * values[6]
					- values[1] * values[3] * values[8]
					- values[0] * values[5] * values[7]);
}

inline void inverse(const DenseMatrix &m, DenseMatrix &inv, double det)
{
	const double *values = m.values();
	inv.resize(m.rows(), m.columns());
	double *invj = inv.values();
	double detJx = 1 / det;
	invj[0] = detJx * (values[8] * values[4] - values[7] * values[5]);
	invj[1] = detJx * (-values[8] * values[1] + values[7] * values[2]);
	invj[2] = detJx * (values[5] * values[1] - values[4] * values[2]);
	invj[3] = detJx * (-values[8] * values[3] + values[6] * values[5]);
	invj[4] = detJx * (values[8] * values[0] - values[6] * values[2]);
	invj[5] = detJx * (-values[5] * values[0] + values[3] * values[2]);
	invj[6] = detJx * (values[7] * values[3] - values[6] * values[4]);
	invj[7] = detJx * (-values[7] * values[0] + values[6] * values[1]);
	invj[8] = detJx * (values[4] * values[0] - values[3] * values[1]);
}

// B =
// dX   0   0
//  0  dY   0
//  0   0  dZ
// dY  dX   0
//  0  dZ  dY
// dZ   0  dX
inline void distribute(DenseMatrix &B, DenseMatrix &dND)
{
	// TODO: block ordering inside B
	eslocal columns = dND.rows() * dND.columns();
	const double *dNDx = dND.values();
	const double *dNDy = dND.values() + dND.columns();
	const double *dNDz = dND.values() + 2 * dND.columns();

	double *v = B.values();

	memcpy(&v[0], dNDx,                               sizeof(double) * dND.columns());
	memcpy(&v[3 * columns + dND.columns()],     dNDx, sizeof(double) * dND.columns());
	memcpy(&v[5 * columns + 2 * dND.columns()], dNDx, sizeof(double) * dND.columns());

	memcpy(&v[1 * columns + dND.columns()],     dNDy, sizeof(double) * dND.columns());
	memcpy(&v[3 * columns],                     dNDy, sizeof(double) * dND.columns());
	memcpy(&v[4 * columns + 2 * dND.columns()], dNDy, sizeof(double) * dND.columns());

	memcpy(&v[2 * columns + 2 * dND.columns()], dNDz, sizeof(double) * dND.columns());
	memcpy(&v[4 * columns + dND.columns()],     dNDz, sizeof(double) * dND.columns());
	memcpy(&v[5 * columns],                     dNDz, sizeof(double) * dND.columns());
}


void Mesh::_assembleElesticity(
		const Element *e,
		size_t part,
		DenseMatrix &Ke,
		DenseMatrix &Me,
		std::vector<double> &fe,
		std::vector<double> &inertia,
		DenseMatrix &C,
		bool dynamic) const
{
	const std::vector<DenseMatrix> &dN = e->dN();
	const std::vector<DenseMatrix> &N = e->N();
	const std::vector<double> &weighFactor = e->weighFactor();

	DenseMatrix coordinates(e->size(), Point::size());
	for (size_t i = 0; i < e->size(); i++) {
		coordinates.values() + i * Point::size() << _coordinates.get(e->node(i), part);
	}

	eslocal Ksize = Point::size() * e->size();
	double detJ;
	DenseMatrix J, invJ, dND, B(C.rows(), Ksize);

	Ke.resize(Ksize, Ksize);
	Ke = 0;
	fe.resize(Ksize);
	fill(fe.begin(), fe.end(), 0);
	if (dynamic) {
		Me.resize(e->size(), e->size());
		Me = 0;
	}

	for (eslocal gp = 0; gp < e->gpSize(); gp++) {
		J.multiply(dN[gp], coordinates);
		detJ = determinant3x3(J);
		inverse(J, invJ, detJ);

		dND.multiply(invJ, dN[gp]);
		distribute(B, dND);

		Ke.multiply(B, C * B, detJ * weighFactor[gp], 1, true);

		for (eslocal i = 0; i < Ksize; i++) {
			fe[i] += detJ * weighFactor[gp] * N[gp](0, i % e->size()) * inertia[i / e->size()];
		}

		if (dynamic) {
			// Me = Me + WF * (DENS * dJ) * (N' * N);
			Me.multiply(N[gp], N[gp], 7.85e-9 * detJ * weighFactor[gp], 1, true);
		}
	}
}

void Mesh::_heat(SparseVVPMatrix<eslocal> &K, SparseVVPMatrix<eslocal> &M, std::vector<double> &f, eslocal part, bool dynamic) const
{
//TODO nK various DOFs per nodes (may be)
	eslocal nK = _coordinates.localSize(part);
	K.resize(nK, nK);
	if (dynamic) {
		M.resize(nK, nK);
	}
	f.resize(nK);

	DenseMatrix Ke, Me;
	std::vector<double> fe;

//	std::vector<double> inertia(3, 0.0);
//	inertia[2] = 9810.0 * 7.85e-9;
//	double ex = 2.1e5;
//	double mi = 0.3;
//	double E = ex / ((1 + mi) * (1 - 2 * mi));
//	DenseMatrix C(6, 6);

//	C(0, 1) = C(0, 2) = C(1, 0) = C(1, 2) = C(2, 0) = C(2, 1) = E * mi;
//	C(0, 0) = C(1, 1) = C(2, 2) = E * (1.0 - mi);
//	C(3, 3) = C(4, 4) = C(5, 5) = E * (0.5 - mi);
	DenseMatrix C(3, 3);
  double kxx = 1.0;
  C(0,0) = kxx;
  C(1,1) = kxx;
  C(2,2) = kxx;

	for (eslocal i = _partPtrs[part]; i < _partPtrs[part + 1]; i++) {
		_assembleHeat(_elements[i], part, Ke, Me, fe, C, dynamic);
		_integrateHeat(_elements[i], K, M, f, Ke, Me, fe, dynamic);
	}



}

void Mesh::_assembleHeat(
		const Element *e,
		size_t part,
		DenseMatrix &Ke,
		DenseMatrix &Me,
		std::vector<double> &fe,
		DenseMatrix &C,
		bool dynamic) const
{
	const std::vector<DenseMatrix> &dN = e->dN();
	const std::vector<DenseMatrix> &N = e->N();
	const std::vector<double> &weighFactor = e->weighFactor();

	DenseMatrix coordinates(e->size(), Point::size());
	for (size_t i = 0; i < e->size(); i++) {
		coordinates.values() + i * Point::size() << _coordinates.get(e->node(i), part);
	}

  double CP = 1.0, rho = 7.85e-9;
	eslocal Ksize =  e->size();
	double detJ;
	DenseMatrix J, invJ, dND;//, B(C.rows(), Ksize);

	Ke.resize(Ksize, Ksize);
	Ke = 0;
	fe.resize(Ksize);
	fill(fe.begin(), fe.end(), 0);
	if (dynamic) {
		Me.resize(e->size(), e->size());
		Me = 0;
	}

	for (eslocal gp = 0; gp < e->gpSize(); gp++) {
		J.multiply(dN[gp], coordinates);
		detJ = determinant3x3(J);
		inverse(J, invJ, detJ);

		dND.multiply(invJ, dN[gp]);
		//distribute(B, dND);
    DenseMatrix &B=dND;


		Ke.multiply(B, C * B, detJ * weighFactor[gp], 1, true);

//		for (eslocal i = 0; i < Ksize; i++) {
//			fe[i] += detJ * weighFactor[gp] * N[gp](0, i % e->size()) * inertia[i / e->size()];
//		}

		if (dynamic) {
			// Me = Me + WF * (DENS * dJ) * (N' * N);
			Me.multiply(N[gp], N[gp], rho * detJ * weighFactor[gp] * CP, 1, true);
		}
	}
}

void Mesh::_integrateHeat(
		const Element *e,
		SparseVVPMatrix<eslocal> &K,
		SparseVVPMatrix<eslocal> &M,
		std::vector<double> &f,
		const DenseMatrix &Ke,
		const DenseMatrix &Me,
		const std::vector<double> &fe,
		bool dynamic
	) const
{
  //
	size_t row, column;
	size_t s = 1;

	for (size_t i = 0; i < s * e->size(); i++) {
		row = s * (e->node(i % e->size())) + i / e->size();
		for (size_t j = 0; j < s * e->size(); j++) {
			column = s * (e->node(j % e->size())) + j / e->size();
			K(row, column) = Ke(i, j);
		}
		//f[row] += fe[i];
	}
	if (!dynamic) {
		return;
	}
	for (size_t i = 0; i < e->size(); i++) {
		row = s * (e->node(i));
		for (size_t j = 0; j < e->size(); j++) {
			column = s * (e->node(j)); //i
			for (size_t k = 0; k < s; k++) {
				M(row + k, column + k) += Me(i, j);
			}
		}
	}
}


void Mesh::_integrateElasticity(
		const Element *e,
		SparseVVPMatrix<eslocal> &K,
		SparseVVPMatrix<eslocal> &M,
		std::vector<double> &f,
		const DenseMatrix &Ke,
		const DenseMatrix &Me,
		const std::vector<double> &fe,
		bool dynamic) const
{
	// Element ordering: xxxx, yyyy, zzzz,...
	// Global ordering:  xyz, xyz, xyz, xyz, ...
	size_t row, column;
	size_t s = Point::size();

	for (size_t i = 0; i < s * e->size(); i++) {
		row = s * (e->node(i % e->size())) + i / e->size();
		for (size_t j = 0; j < s * e->size(); j++) {
			column = s * (e->node(j % e->size())) + j / e->size();
			K(row, column) = Ke(i, j);
		}
		f[row] += fe[i];
	}
	if (!dynamic) {
		return;
	}
	for (size_t i = 0; i < e->size(); i++) {
		row = s * (e->node(i));
		for (size_t j = 0; j < e->size(); j++) {
			column = s * (e->node(j)); //i
			for (size_t k = 0; k < s; k++) {
				M(row + k, column + k) += Me(i, j);
			}
		}
	}
}

void Mesh::partitiate(eslocal parts, eslocal fixPoints)
{
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
	}
}

void Mesh::computeFixPoints(eslocal fixPoints)
{
	_fixPoints.resize(parts() * fixPoints);

#ifndef DEBUG
	cilk_for (eslocal i = 0; i < parts(); i++) {
#else
	for (eslocal i = 0; i < parts(); i++) {
#endif
		eslocal *eSubPartition = getPartition(_partPtrs[i], _partPtrs[i + 1], fixPoints);

		for (eslocal j = 0; j < fixPoints; j++) {
			_fixPoints[i * fixPoints + j] = getCentralNode(_partPtrs[i], _partPtrs[i + 1], eSubPartition, i, j);
		}

		delete[] eSubPartition;
	}
}

eslocal* Mesh::getPartition(eslocal first, eslocal last, eslocal parts) const
{
	if (parts == 1) {
		eslocal *ePartition = new eslocal[last - first];
		for (eslocal i = first; i < last; i++) {
			ePartition[i] = 0;
		}
		return ePartition;
	}
	// INPUTS
	eslocal ncommon, eSize, nSize, *e, *n, options[METIS_NOPTIONS];

	// OUTPUTS
	eslocal objval, *ePartition, *nPartition;

	// FILL INPUT VARIABLES
	////////////////////////////////////////////////////////////////////////////

	// number of common nodes to be neighbor
	ncommon = 2;

	// There is probably BUG in METIS numbering or I do not understand documentation.
	// The solution is increase the size of 'nodesCount' and keep the default numbering
	//options[METIS_OPTION_NUMBERING] = coordinates.getNumbering();
	METIS_SetDefaultOptions(options);

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
	for (eslocal i = first, index = 0; i < last; i++, index++) {
		_elements[i]->fillNodes(n + e[index]);
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
		_elements[e]->setLocalIndices(nodeMap);
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

void Mesh::saveBasis(
		std::ofstream &vtk,
		std::vector<std::vector<eslocal> > &l2g_vec,
		const Boundaries &lBoundaries,
		const Boundaries &gBoundaries,
		double subDomainShrinking,
		double clusterShrinking) const
{
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

	Point cCenter;
	for (size_t i = 0; i < _coordinates.size(); i++) {
		cCenter += _coordinates[i];
	}
	cCenter /= _coordinates.size();

	vtk << "POINTS " << n_points << " float\n";
	for (size_t d = 0; d < nSubClst; d++) {
		Point sCenter;
		for (size_t c = 0; c < l2g_vec[d].size(); c++) {
			sCenter += _coordinates[l2g_vec[d][c]];
		}
		sCenter /= l2g_vec[d].size();

		for (size_t i = 0; i < l2g_vec[d].size(); i++) {
			Point x = _coordinates[l2g_vec[d][i]];
			x = sCenter + (x - sCenter) * subDomainShrinking;
			x = cCenter + (x - cCenter) * clusterShrinking;
			vtk << x << "\n";
		}
	}

	size_t size = 0;
	for (size_t i = 0; i < _elements.size(); i++) {
		size += _elements[i]->size() + 1;
	}

	if (getFixPointsCount()) {
		size += parts() + _fixPoints.size();
		vtk << "CELLS " << _elements.size() + parts() << " " << size << "\n";
	} else {
		vtk << "CELLS " << _elements.size() << " " << size << "\n";
	}

	size_t i = 0;
	for (size_t part = 0; part + 1 < _partPtrs.size(); part++) {
		for (eslocal ii = 0; ii < _partPtrs[part + 1] - _partPtrs[part]; ii++) {
			vtk << _elements[i]->size();
			for (size_t j = 0; j < _elements[i]->size(); j++) {
				vtk << " " << _elements[i]->node(j) + cnt;
			}
			vtk << "\n";
			i++;
		}
		cnt += l2g_vec[part].size();
	}
	if (getFixPointsCount()) {
		cnt = 0;
		for (size_t part = 0; part < parts(); part++) {
			vtk << getFixPointsCount();
			for (size_t i = 0; i < getFixPointsCount(); i++) {
				vtk << " " << _fixPoints[part * getFixPointsCount() + i] + cnt;
			}
			vtk << "\n";
			cnt += l2g_vec[part].size();
		}
	}

	vtk << "\n";
	if (getFixPointsCount()) {
		vtk << "CELL_TYPES " << _elements.size() + parts() << "\n";
	} else {
		vtk << "CELL_TYPES " << _elements.size() << "\n";
	}
	for (size_t i = 0; i < _elements.size(); i++) {
		vtk << _elements[i]->vtkCode() << "\n";
	}
	if (getFixPointsCount()) {
		for (size_t p = 0; p + 1 < _partPtrs.size(); p++) {
			vtk << "2\n";
		}
	}

	vtk << "\n";
	if (getFixPointsCount()) {
		vtk << "CELL_DATA " << _elements.size() + parts() << "\n";
	} else {
		vtk << "CELL_DATA " << _elements.size() << "\n";
	}
	vtk << "SCALARS decomposition int 1\n";
	vtk << "LOOKUP_TABLE decomposition\n";
	for (size_t part = 0; part + 1 < _partPtrs.size(); part++) {
		for (eslocal i = 0; i < _partPtrs[part + 1] - _partPtrs[part]; i++) {
			vtk << part << "\n";
		}
	}
	if (getFixPointsCount()) {
		for (size_t part = 0; part + 1 < _partPtrs.size(); part++) {
			vtk << parts() << "\n";
		}
	}
}

void Mesh::saveVTK(
		const char* filename,
		std::vector<std::vector<double> > &displacement,
		std::vector<std::vector<eslocal> > &l2g_vec,
		const Boundaries &lBoundaries,
		const Boundaries &gBoundaries,
		double subDomainShrinking,
		double clusterShrinking) const
{
	std::ofstream vtk(filename, std::ios::out | std::ios::trunc);
	saveBasis(vtk, l2g_vec, lBoundaries, gBoundaries, subDomainShrinking, clusterShrinking);

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

void Mesh::saveVTK(
		const char* filename,
		const Boundaries &localBoundaries,
		double subDomainShrinking,
		double clusterShrinking) const
{
	std::ofstream vtk;

	size_t cSize = 0;
	for (size_t p = 0; p < parts(); p++) {
		cSize += _coordinates.localToCluster(p).size();
	}

	vtk.open(filename, std::ios::out | std::ios::trunc);
	vtk << "# vtk DataFile Version 3.0\n";
	vtk << "Test\n";
	vtk << "ASCII\n\n";
	vtk << "DATASET UNSTRUCTURED_GRID\n";
	vtk << "POINTS " << cSize << " float\n";

	Point cCenter;
	for (size_t i = 0; i < _coordinates.size(); i++) {
		cCenter += _coordinates[i];
	}
	cCenter /= _coordinates.size();

	for (size_t p = 0; p < parts(); p++) {

		Point sCenter;
		for (size_t i = 0; i < _coordinates.localSize(p); i++) {
			sCenter += _coordinates.get(i, p);
		}
		sCenter /= _coordinates.localSize(p);

		for (size_t i = 0; i < _coordinates.localSize(p); i++) {
			Point x = _coordinates.get(i, p);
			x = sCenter + (x - sCenter) * subDomainShrinking;
			x = cCenter + (x - cCenter) * clusterShrinking;
			vtk << x << "\n";
		}
	}
	vtk << "\n";

	size_t size = 0;
	for (size_t i = 0; i < _elements.size(); i++) {
		size += _elements[i]->size() + 1;
	}

	size_t corners = 0;
	for (size_t i = 0; i < localBoundaries.size(); i++) {
		if (localBoundaries.isCorner(i)) {
			corners += localBoundaries[i].size();
		}
	}

	vtk << "CELLS " << _elements.size() + 1 << " " << size + corners  + 1 << "\n";

	size_t offset = 0;
	for (size_t p = 0; p < parts(); p++) {
		for (size_t i = _partPtrs[p]; i < _partPtrs[p + 1]; i++) {
			vtk << _elements[i]->size();
			for (size_t j = 0; j < _elements[i]->size(); j++) {
				vtk << " " << _elements[i]->node(j) + offset;
			}
			vtk << "\n";
		}
		offset += _coordinates.localToCluster(p).size();
	}

	offset = 0;
	vtk << corners;
	for (size_t p = 0; p < parts(); p++) {
		for (size_t i = 0; i < localBoundaries.size(); i++) {
			if (localBoundaries.isCorner(i) && localBoundaries[i].count(p)) {
				vtk << " " << _coordinates.localIndex(i, p) + offset;
			}
		}
		offset += _coordinates.localSize(p);
	}

	vtk << "\n";
	vtk << "CELL_TYPES " << _elements.size() + 1 << "\n";

	for (size_t i = 0; i < _elements.size(); i++) {
		vtk << _elements[i]->vtkCode() << "\n";
	}
	vtk << "2\n";

	vtk << "\n";
	vtk << "CELL_DATA " << _elements.size() + 1 << "\n";

	vtk << "SCALARS decomposition int 1\n";
	vtk << "LOOKUP_TABLE decomposition\n";
	for (size_t part = 0; part + 1 < _partPtrs.size(); part++) {
		for (eslocal i = 0; i < _partPtrs[part + 1] - _partPtrs[part]; i++) {
			vtk << part << "\n";
		}
	}
	vtk << parts() << "\n";

	vtk.close();
}

void Mesh::saveVTK(
		const char* filename,
		double subDomainShrinking,
		double clusterShrinking) const
{
	std::ofstream vtk;

	size_t cSize = 0;
	for (size_t p = 0; p < parts(); p++) {
		cSize += _coordinates.localToCluster(p).size();
	}

	vtk.open(filename, std::ios::out | std::ios::trunc);
	vtk << "# vtk DataFile Version 3.0\n";
	vtk << "Test\n";
	vtk << "ASCII\n\n";
	vtk << "DATASET UNSTRUCTURED_GRID\n";
	vtk << "POINTS " << cSize << " float\n";

	Point cCenter;
	for (size_t i = 0; i < _coordinates.size(); i++) {
		cCenter += _coordinates[i];
	}
	cCenter /= _coordinates.size();


	for (size_t p = 0; p + 1 < _partPtrs.size(); p++) {

		Point sCenter;
		for (size_t i = 0; i < _coordinates.localSize(p); i++) {
			sCenter += _coordinates.get(i, p);
		}
		sCenter /= _coordinates.localSize(p);

		for (size_t i = 0; i < _coordinates.localSize(p); i++) {
			Point x = _coordinates.get(i, p);
			x = sCenter + (x - sCenter) * subDomainShrinking;
			x = cCenter + (x - cCenter) * clusterShrinking;
			vtk << x << "\n";
			// TODO: Catalyst
			// Grid::Points x contains all three coordinates (x, y and z)
			// std::vector<double> Points;
			// pointArray->SetArray(grid.GetPointsArray(), static_cast<vtkIdType>(grid.GetNumberOfPoints()*3), 1); in FEAdaptor.cxx line 28
		}
	}
	vtk << "\n";

	size_t size = 0;
	for (size_t i = 0; i < _elements.size(); i++) {
		size += _elements[i]->size() + 1;
	}

	if (getFixPointsCount()) {
		size += parts() + _fixPoints.size();
		vtk << "CELLS " << _elements.size() + parts() << " " << size << "\n";
	} else {
		vtk << "CELLS " << _elements.size() << " " << size << "\n";
	}
	size_t offset = 0;
	for (size_t p = 0; p + 1 < _partPtrs.size(); p++) {
	// TODO: Catalyst : create the cells
	// create the cells
	// size_t numCells = _elements->size();
	// VTKGrid->Allocate(static_cast<vtkIdType>(size));
		for (size_t i = _partPtrs[p]; i < _partPtrs[p + 1]; i++) {
			vtk << _elements[i]->size();
			//vtkIdType tmp[_elements[i]->size()];
			for (size_t j = 0; j < _elements[i]->size(); j++) {
				vtk << " " << _elements[i]->node(j) + offset;
				//vtkIdType tmp[j] = _elements[i]->node(j) + offset;
			}
			//VTKGrid->InsertNextCell(_elements[i]->vtkCode(), 8, tmp);
			vtk << "\n";
		}
		offset += _coordinates.localToCluster(p).size();
	}

	if (getFixPointsCount()) {
		offset = 0;
		for (size_t p = 0; p < parts(); p++) {
			vtk << getFixPointsCount();
			for (size_t i = 0; i < getFixPointsCount(); i++) {
				vtk << " " << _fixPoints[p * getFixPointsCount() + i] + offset;
			}
			vtk << "\n";
			offset += _coordinates.localToCluster(p).size();
		}
	}

	vtk << "\n";
	if (getFixPointsCount()) {
		vtk << "CELL_TYPES " << _elements.size() + parts() << "\n";
	} else {
		vtk << "CELL_TYPES " << _elements.size() << "\n";
	}
	for (size_t i = 0; i < _elements.size(); i++) {
		vtk << _elements[i]->vtkCode() << "\n";
	}
	if (getFixPointsCount()) {
		for (size_t p = 0; p + 1 < _partPtrs.size(); p++) {
			vtk << "2\n";
		}
	}

	vtk << "\n";
	if (getFixPointsCount()) {
		vtk << "CELL_DATA " << _elements.size() + parts() << "\n";
	} else {
		vtk << "CELL_DATA " << _elements.size() << "\n";
	}

	vtk << "SCALARS decomposition int 1\n";
	vtk << "LOOKUP_TABLE decomposition\n";
	for (size_t part = 0; part + 1 < _partPtrs.size(); part++) {
		for (eslocal i = 0; i < _partPtrs[part + 1] - _partPtrs[part]; i++) {
			vtk << part << "\n";
		}
	}
	if (getFixPointsCount()) {
		for (size_t part = 0; part + 1 < _partPtrs.size(); part++) {
			vtk << parts() << "\n";
		}
	}
	vtk << "\n";

	vtk.close();
}

bool isOuterFace(
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

eslocal isOnBoundary(
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
	surface._partPtrs.push_back(surface._elements.size());
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
		surface.computeLocalIndices(surface._partPtrs.size() - 1);
	}
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

void Mesh::computeCorners(Boundaries &boundaries, eslocal number, bool corners, bool edges, bool faces) const
{
	if (parts() < 1) {
		std::cerr << "Internal error: _partPtrs.size()\n";
		exit(EXIT_FAILURE);
	}
	if (!corners && !edges && !faces) {
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

#ifndef DEBUG
	cilk_for (eslocal i = 0; i < parts(); i++) {
#else
	for (eslocal i = 0; i < parts(); i++) {
#endif
		// compute number of elements and fill used nodes
		for (eslocal j = _partPtrs[i]; j < _partPtrs[i + 1]; j++) {
			for (size_t k = 0; k < _elements[j]->faces(); k++) {
				std::vector<eslocal> face = _elements[j]->getFace(k);
				for (size_t f = 0; f < face.size(); f++) {
					face[f] = _coordinates.clusterIndex(face[f], i);
				}
				if (isOnBoundary(nodesElements, face, _partPtrs, 1) == j) {
					commonFaces.push_back(face);
				}
			}
		}
	}

	// create mesh from common faces
	Mesh cfm;
	Element *e;
	if (faces) {
		cfm.coordinates() = _coordinates;
		cfm._elements.reserve(commonFaces.size());
	}
	std::vector<std::vector<eslocal> > nodesFaces(_coordinates.size());
	for (size_t i = 0; i < commonFaces.size(); i++) {
		if (commonFaces[i].size() == 3) {
			e = new Triangle(&commonFaces[i][0]);
		}
		if (commonFaces[i].size() == 4) {
			e = new Square(&commonFaces[i][0]);
		}
		if (faces) {
			cfm._elements.push_back(e);
		}
		// preparation for corners on edges
		for (size_t j = 0; j < e->size(); j++) {
			std::vector<eslocal> neigh = e->getNeighbours(j);
			eslocal index = e->node(j);
			for (size_t n = 0; n < neigh.size(); n++) {
				if (index < neigh[n]) {
					neighbours[index].insert(neigh[n]);
				}
			}
			nodesFaces[index].push_back(i);
		}
	}
	if (faces) {
		cfm._partPtrs.back() = cfm._elements.size();
		cfm.computeLocalIndices(cfm._partPtrs.size() - 1);
		cfm.computeFixPoints(number);
		for (size_t i = 0; i < cfm.getFixPointsCount(); i++) {
			boundaries.setCorner(cfm.coordinates().clusterIndex(cfm.getFixPoints()[i], 0));
		}
	}

	Mesh clm;
	if (edges) {
		clm.coordinates() = _coordinates;
	}

	std::vector<eslocal> result;
	std::vector<eslocal>::iterator end;
	std::vector<eslocal> pairs;
	std::vector<eslocal> pair(2);
	std::vector<eslocal> nodeCounter(_coordinates.size(), 0);

	std::set<eslocal>::iterator it;
	for (size_t i = 0; i < neighbours.size(); i++) {
		for (it = neighbours[i].begin(); it != neighbours[i].end(); ++it) {
			pair[0] = i;
			pair[1] = *it;
			result.resize(boundaries[i].size());
			end = std::set_intersection(
				boundaries[i].begin(), boundaries[i].end(),
				boundaries[*it].begin(), boundaries[*it].end(),
				result.begin());
			if (end - result.begin() >= 3 && isOnBoundary(nodesElements, pair, _partPtrs, 2)) {
				pairs.push_back(pair[0]);
				pairs.push_back(pair[1]);
				nodeCounter[pair[0]]++;
				nodeCounter[pair[1]]++;
				continue;
			}
			if (isOuterFace(nodesFaces, pair)) {
				pairs.push_back(pair[0]);
				pairs.push_back(pair[1]);
				nodeCounter[pair[0]]++;
				nodeCounter[pair[1]]++;
			}
		}
	}

	if (edges) {
		clm._elements.reserve(pairs.size() / 2);
	}

	for (size_t j = 0; j < pairs.size(); j += 2) {
		if (corners) {
			if (nodeCounter[pairs[j]] > 3 || nodeCounter[pairs[j + 1]] > 3) {
				if (nodeCounter[pairs[j]] > 3) {
					boundaries.setCorner(pairs[j]);
				}
				if (nodeCounter[pairs[j + 1]] > 3) {
					boundaries.setCorner(pairs[j + 1]);
				}
				continue;
			}
		}
		if (edges) {
			clm._elements.push_back(new Line(&pairs[j]));
		}
	}
	if (edges) {
		clm._partPtrs.back() = clm._elements.size();
		clm.computeLocalIndices(clm._partPtrs.size() - 1);
		clm.computeFixPoints(number);
		for (size_t i = 0; i < clm.getFixPointsCount(); i++) {
			boundaries.setCorner(clm.coordinates().clusterIndex(clm.getFixPoints()[i], 0));
		}
	}
}

void SurfaceMesh::elasticity(DenseMatrix &K, size_t part) const
{
	eslocal nK = Point::size() * _coordinates.localSize(part);
	eslocal eSize = _partPtrs[part + 1] - _partPtrs[part];
	K.resize(nK, nK);
	std::vector<double> nodes(nK);
	std::vector<eslocal> elems(3 * eSize);

	for (size_t i = 0; i < _coordinates.localSize(part); i++) {
		&nodes[i * Point::size()] << _coordinates.get(i, part);
	}
	for (size_t i = _partPtrs[part], index = 0; i < _partPtrs[part + 1];
			i++, index++) {
		// TODO: various data types int32_t and int64_t
		// _elements[i]->fillNodes(&elems[3 * i]); CANNOT be used
		for (size_t j = 0; j < _elements[i]->size(); j++) {
			elems[3 * index + j] = _elements[i]->node(j);
		}
	}

	bem4i::getLameSteklovPoincare(
			K.values(),
			_coordinates.localSize(part),
			&nodes[0],
			eSize,
			&elems[0],
			0.33,			// nu
			1.0e5,			// E
			3,				// order near
			4,				// order far
			false			// verbose
			);
}

void SurfaceMesh::integrateUpperFaces(std::vector<double> &f, size_t part)
{
	double hight_z = 29.99999999;
	Point p0, p1, p2, v10, v20;
	double Area_h;

	for (size_t i = _partPtrs[part]; i < _partPtrs[part + 1]; i++) {
		bool flag_edgeOnTop = true;
		for (size_t j = 0; j < _elements[i]->size(); j++) {
			if (_coordinates.get(_elements[i]->node(j), part).z < hight_z) {
				flag_edgeOnTop = false;
				continue;
			}
		}
		if (flag_edgeOnTop) {
			p0 = _coordinates.get(_elements[i]->node(0), part);
			p1 = _coordinates.get(_elements[i]->node(1), part);
			p2 = _coordinates.get(_elements[i]->node(2), part);
			v10 = p1 - p0;
			v20 = p2 - p0;

			Area_h = 0.5
					* (v10.y * v20.z - v20.y * v10.z + v20.x * v10.z
					-  v10.x * v20.z + v10.x * v20.y - v20.x * v10.y);

			for (size_t k = 0; k < 3; k++) {
				f[3 * _elements[i]->node(k) + 2] += (1. / 3.) * Area_h;
			}
		}
	}
}


void Mesh::saveData()
{
	eslocal value;
	esglobal index;

	for (size_t p = 0; p < parts(); p++) {
		std::stringstream ss;
		ss << "mesh" << p << ".dat";
		std::ofstream os(ss.str().c_str(), std::ofstream::binary | std::ofstream::trunc);

		// save elements
		value = _partPtrs[p + 1] - _partPtrs[p];
		os.write(reinterpret_cast<const char*>(&value), sizeof(eslocal));
		for (eslocal e = _partPtrs[p]; e < _partPtrs[p + 1]; e++) {
			os << *(_elements[e]);
		}

		// save coordinates
		value = _coordinates.localSize(p);
		os.write(reinterpret_cast<const char*>(&value), sizeof(eslocal));
		const std::vector<eslocal> &l2c = _coordinates.localToCluster(p);
		for (eslocal i = 0; i < l2c.size(); i++) {
			index = _coordinates.globalIndex(i, p);
			os.write(reinterpret_cast<const char*>(&index), sizeof(esglobal));
			const Point &point = _coordinates.get(i, p);
			os.write(reinterpret_cast<const char*>(&point), Point::size() * sizeof(double));
		}

		// save coordinates' properties
		for (size_t i = 0; i < _coordinates.propertiesSize(); i++) {
			const std::map<eslocal, double> &property = _coordinates.property(i).values();
			eslocal counter = 0;
			const std::vector<eslocal> &l2c = _coordinates.localToCluster(p);
			for (size_t j = 0; j < l2c.size(); j++) {
				if (property.find(l2c[j]) != property.end()) {
					counter++;
				}
			}
			os.write(reinterpret_cast<const char*>(&counter), sizeof(eslocal));
			for (eslocal j = 0; j < l2c.size(); j++) {
				if (property.find(l2c[j]) != property.end()) {
					os.write(reinterpret_cast<const char*>(&j), sizeof(eslocal));
					os.write(reinterpret_cast<const char*>(&property.find(l2c[j])->second), sizeof(double));
				}
			}
		}
	}
}

void Mesh::loadData(const char *filename)
{
	std::ifstream is(filename, std::ifstream::binary);

	// reset parameters
	eslocal size, type, cIndex;
	esglobal index;
	double value;
	Point point;

	is.read(reinterpret_cast<char *>(&size), sizeof(eslocal));
	_elements.clear();
	_partPtrs.resize(2);
	_partPtrs[1] = size;
	_fixPoints.clear();
	_elements.reserve(size);
	_coordinates.clear();

	// load elements
	for (eslocal i = 0; i < size; i++) {
		is.read(reinterpret_cast<char *>(&type), sizeof(eslocal));
		switch(type) {
		case Tetrahedron4VTKCode:
			_elements.push_back(new Tetrahedron4(is));
			break;
		case Tetrahedron10VTKCode:
			_elements.push_back(new Tetrahedron10(is));
			break;
		case Pyramid5VTKCode:
			_elements.push_back(new Pyramid5(is));
			break;
		case Pyramid13VTKCode:
			_elements.push_back(new Pyramid13(is));
			break;
		case Prisma6VTKCode:
			_elements.push_back(new Prisma6(is));
			break;
		case Prisma15VTKCode:
			_elements.push_back(new Prisma15(is));
			break;
		case Hexahedron8VTKCode:
			_elements.push_back(new Hexahedron8(is));
			break;
		case Hexahedron20VTKCode:
			_elements.push_back(new Hexahedron20(is));
			break;
		}
	}

	// load coordinates
	is.read(reinterpret_cast<char *>(&size), sizeof(eslocal));
	_coordinates.reserve(size);
	for (eslocal i = 0; i < size; i++) {
		is.read(reinterpret_cast<char *>(&index), sizeof(esglobal));
		is.read(reinterpret_cast<char *>(&point), Point::size() * sizeof(double));
		_coordinates.add(point, i, index);
	}

	// load coordinates' properties
	for (size_t i = 0; i < _coordinates.propertiesSize(); i++) {
		CoordinatesProperty &property = _coordinates.property(i);
		is.read(reinterpret_cast<char *>(&size), sizeof(eslocal));
		for (eslocal j = 0; j < size; j++) {
			is.read(reinterpret_cast<char *>(&cIndex), sizeof(eslocal));
			is.read(reinterpret_cast<char *>(&value), sizeof(double));
			property[cIndex] = value;
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

