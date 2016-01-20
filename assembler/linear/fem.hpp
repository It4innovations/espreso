
#include "linear.h"

namespace assembler {

static double determinant3x3(DenseMatrix &m)
{
	const double *values = m.values();
	return fabs(
		values[0] * values[4] * values[8] +
		values[1] * values[5] * values[6] +
		values[2] * values[3] * values[7] -
		values[2] * values[4] * values[6] -
		values[1] * values[3] * values[8] -
		values[0] * values[5] * values[7]
   );
}

static void inverse(const DenseMatrix &m, DenseMatrix &inv, double det)
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
static void distribute(DenseMatrix &B, DenseMatrix &dND)
{
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


template <>
void Linear<FEM>::KeMefe(
		DenseMatrix &Ke, DenseMatrix &Me, std::vector<double> &fe,
		DenseMatrix &Ce, const mesh::Element *e, size_t part, bool dynamics)
{
	const std::vector<DenseMatrix> &dN = e->dN();
	const std::vector<DenseMatrix> &N = e->N();
	const std::vector<double> &weighFactor = e->weighFactor();
	std::vector<double> inertia;
	this->inertia(inertia);

	DenseMatrix coordinates(e->size(), mesh::Point::size());
	for (size_t i = 0; i < e->size(); i++) {
		coordinates.values() + i * mesh::Point::size() << _input.mesh.coordinates().get(e->node(i), part);
	}

	eslocal Ksize = e->size() * this->DOFs();
	double detJ;
	DenseMatrix J, invJ, dND;

	Ke.resize(Ksize, Ksize);
	Ke = 0;
	fe.resize(Ksize);
	fill(fe.begin(), fe.end(), 0);
	if (dynamics) {
		Me.resize(e->size(), e->size());
		Me = 0;
	}

	for (eslocal gp = 0; gp < e->gpSize(); gp++) {
		J.multiply(dN[gp], coordinates);
		detJ = determinant3x3(J);
		inverse(J, invJ, detJ);

		dND.multiply(invJ, dN[gp]);

		// TODO: make it more general
		if (this->DOFs() == 3) {
			DenseMatrix B(Ce.rows(), Ksize);
			distribute(B, dND);
			Ke.multiply(B, Ce * B, detJ * weighFactor[gp], 1, true);
		} else {
			Ke.multiply(dND, Ce * dND, detJ * weighFactor[gp], 1, true);
		}

		for (eslocal i = 0; i < Ksize; i++) {
			fe[i] += detJ * weighFactor[gp] * N[gp](0, i % e->size()) * inertia[i / e->size()];
		}

		if (dynamics) {
			// Me = Me + WF * (DENS * dJ) * (N' * N);
			Me.multiply(N[gp], N[gp], this->rho() * detJ * weighFactor[gp] * this->CP(), 1, true);
		}
	}
}

template <>
void Linear<FEM>::integrate(
			DenseMatrix &Ke, DenseMatrix &Me, std::vector<double> &fe,
			SparseVVPMatrix<eslocal> &K, SparseVVPMatrix<eslocal> &M, std::vector<double> &f,
			const mesh::Element *e, bool dynamics)
{
	// Element ordering: xxxx, yyyy, zzzz,...
	// Global ordering:  xyz, xyz, xyz, xyz, ...
	size_t row, column;
	size_t s = this->DOFs();

	for (size_t i = 0; i < s * e->size(); i++) {
		row = s * (e->node(i % e->size())) + i / e->size();
		for (size_t j = 0; j < s * e->size(); j++) {
			column = s * (e->node(j % e->size())) + j / e->size();
			K(row, column) = Ke(i, j);
		}
		f[row] += fe[i];
	}
	if (!dynamics) {
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


template <>
void Linear<FEM>::KMf(size_t part, bool dynamics)
{
	SparseVVPMatrix<eslocal> _K;
	SparseVVPMatrix<eslocal> _M;
	eslocal nK = _input.mesh.coordinates().localSize(part) * this->DOFs();
	_K.resize(nK, nK);
	if (dynamics) {
		_M.resize(nK, nK);
	}
	_f[part].resize(nK);

	DenseMatrix Ke, Me, Ce;
	std::vector<double> fe;

	this->C(Ce);

	const std::vector<eslocal> &partition = _input.mesh.getPartition();
	const std::vector<mesh::Element*> &elements = _input.mesh.getElements();
	for (eslocal i = partition[part]; i < partition[part + 1]; i++) {
		KeMefe(Ke, Me, fe, Ce, elements[i], part, dynamics);
		integrate(Ke, Me, fe, _K, _M, _f[part], elements[i], dynamics);
	}

	// TODO: make it direct
	SparseCSRMatrix<eslocal> csrK = _K;
	SparseCSRMatrix<eslocal> csrM = _M;

	this->_K[part] = csrK;

	std::stringstream ss;
	ss << "K" << part << ".txt";
	std::ofstream os(ss.str().c_str());
	os << csrK;
	os.close();

	this->_M[part] = csrM;
}

template <>
void Linear<FEM>::T(size_t part)
{
	eslocal nT = _input.mesh.coordinates().localSize(part) * this->DOFs();
	SparseDOKMatrix<eslocal> _T(nT, nT);
	for (size_t i = 0; i < nT; i++) {
		_T(i, i) = 1;
	}

	const mesh::Coordinates& coords = _input.mesh.coordinates();
	const mesh::Boundaries& boundary = _input.mesh.subdomainBoundaries();

	for (size_t i = 0; i < coords.localSize(part); i++) {
		if (boundary.isAveraging(coords.clusterIndex(i, part))) {
			const std::vector<eslocal>& av = boundary.averaging(coords.clusterIndex(i, part));
			for (size_t a = 0; a < av.size(); a++) {
				eslocal j = coords.localIndex(av[a], part);
				for (int d = 0; d < this->DOFs(); d++) {
					_T(i * this->DOFs() + d, j * this->DOFs() + d) = -1;
					_T(j * this->DOFs() + d, i * this->DOFs() + d) = 1;
				}
			}
		}
	}

	SparseCSRMatrix<eslocal> tmpT = _T;

	std::stringstream ss;
	ss << "T" << part << ".txt";
	std::ofstream os(ss.str().c_str());
	os << tmpT;
	os.close();

	this->_T[part] = tmpT;
}

template <>
void Linear<FEM>::initSolver()
{
	_lin_solver.init(
		_input.mesh,
		_K,
		_T,
		_globalB,
		_localB,
		_lambda_map_sub_B1,
		_lambda_map_sub_B0,
		_lambda_map_sub_clst,
		_B1_duplicity,
		_f,
		_vec_c,
		_input.mesh.getFixPoints(),
		_neighClusters
	);
}

template <>
void Linear<FEM>::RHS()
{
	const std::map<eslocal, double> &forces_x = this->_input.mesh.coordinates().property(mesh::FORCES_X).values();
	const std::map<eslocal, double> &forces_y = this->_input.mesh.coordinates().property(mesh::FORCES_Y).values();
	const std::map<eslocal, double> &forces_z = this->_input.mesh.coordinates().property(mesh::FORCES_Z).values();

	for (size_t p = 0; p < this->_input.mesh.parts(); p++) {
		const std::vector<eslocal> &l2g = this->_input.mesh.coordinates().localToCluster(p);
		for (eslocal i = 0; i < l2g.size(); i++) {
			size_t n = this->_input.mesh.subdomainBoundaries()[l2g[i]].size();
			if (forces_x.find(l2g[i]) != forces_x.end()) {
				_f[p][this->DOFs() * i + 0] = forces_x.at(l2g[i]) / n;
			}
			if (forces_y.find(l2g[i]) != forces_y.end()) {
				_f[p][this->DOFs() * i + 1] = forces_y.at(l2g[i]) / n;
			}
			if (forces_z.find(l2g[i]) != forces_z.end()) {
				_f[p][this->DOFs() * i + 2] = forces_z.at(l2g[i]) / n;
			}
		}
		std::stringstream ss;
		ss << "f" << p << ".txt";
		std::ofstream os(ss.str().c_str());
		os << _f[p];
		os.close();
	}

}

}
