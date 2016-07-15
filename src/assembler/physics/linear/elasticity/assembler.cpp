
#include "assembler.h"

using namespace espreso;

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

static void processElement(DenseMatrix &Ke, std::vector<double> &fe, const espreso::Mesh &mesh, size_t subdomain, const Element* element)
{
	DenseMatrix Ce(6, 6), coordinates, J, invJ, dND, B;
	std::vector<double> inertia(3, 0);
	double detJ;

	const Material &material = mesh.materials()[element->getParam(Element::MATERIAL)];
	const std::vector<DenseMatrix> &dN = element->dN();
	const std::vector<DenseMatrix> &N = element->N();
	const std::vector<double> &weighFactor = element->weighFactor();

	// TODO: set the omega from example
	Point omega(50, 50, 0);

	double ex = material.youngModulus;
	double mi = material.poissonRatio;
	double E = ex / ((1 + mi) * (1 - 2 * mi));
	Ce(0, 1) = Ce(0, 2) = Ce(1, 0) = Ce(1, 2) = Ce(2, 0) = Ce(2, 1) = E * mi;
	Ce(0, 0) = Ce(1, 1) = Ce(2, 2) = E * (1.0 - mi);
	Ce(3, 3) = Ce(4, 4) = Ce(5, 5) = E * (0.5 - mi);

	inertia[0] = inertia[1] = 0;
	inertia[2] = 9.8066 * material.density;

	coordinates.resize(element->size(), 3);

	Point mid;
	for (size_t i = 0; i < element->size(); i++) {
		coordinates(i, 0) = mesh.coordinates().get(element->node(i), subdomain).x;
		coordinates(i, 1) = mesh.coordinates().get(element->node(i), subdomain).y;
		coordinates(i, 2) = mesh.coordinates().get(element->node(i), subdomain).z;
		mid += mesh.coordinates().get(element->node(i), subdomain);
	}
	mid /= element->size();

	eslocal Ksize = 3 * element->size();

	Ke.resize(Ksize, Ksize);
	Ke = 0;
	fe.resize(Ksize);
	fill(fe.begin(), fe.end(), 0);

	double rotation[3] = { mid.x * omega.x * omega.x, mid.y * omega.y * omega.y, mid.z * omega.z * omega.z };

	for (eslocal gp = 0; gp < element->gpSize(); gp++) {
		J.multiply(dN[gp], coordinates);
		detJ = determinant3x3(J);
		inverse(J, invJ, detJ);

		dND.multiply(invJ, dN[gp]);
		B.resize(Ce.rows(), Ksize);
		distribute(B, dND);
		Ke.multiply(B, Ce * B, detJ * weighFactor[gp], 1, true);

		for (eslocal i = 0; i < Ksize; i++) {
			// TODO: set rotation from example
			fe[i] += detJ * weighFactor[gp] * N[gp](0, i % element->size()) * inertia[i / element->size()];
			//fe[i] += detJ * weighFactor[gp] * N[gp](0, i % e->size()) * 7850 * rotation[i / e->size()];
		}
	}
}

static void analyticsKernels(SparseMatrix &R1, const Coordinates &coordinates, size_t subdomain)
{
	size_t nodes = coordinates.localSize(subdomain);
	R1.rows = 3 * nodes;
	R1.cols = 6;
	R1.nnz = R1.rows * R1.cols;
	R1.type = 'G';

	R1.dense_values.reserve(R1.nnz);

	for (size_t c = 0; c < 3; c++) {
		std::vector<double> kernel = { 0, 0, 0 };
		kernel[c] = 1;
		for (size_t i = 0; i < nodes; i++) {
			R1.dense_values.insert(R1.dense_values.end(), kernel.begin(), kernel.end());
		}
	}

	for (size_t i = 0; i < coordinates.localSize(subdomain); i++) {
		const Point &p = coordinates.get(i, subdomain);
		R1.dense_values.push_back(-p.y);
		R1.dense_values.push_back( p.x);
		R1.dense_values.push_back(   0);
	}

	for (size_t i = 0; i < coordinates.localSize(subdomain); i++) {
		const Point &p = coordinates.get(i, subdomain);
		R1.dense_values.push_back(-p.z);
		R1.dense_values.push_back(   0);
		R1.dense_values.push_back( p.x);
	}

	for (size_t i = 0; i < coordinates.localSize(subdomain); i++) {
		const Point &p = coordinates.get(i, subdomain);
		R1.dense_values.push_back(   0);
		R1.dense_values.push_back(-p.z);
		R1.dense_values.push_back( p.y);
	}
}

static void analyticsRegMat(SparseMatrix &K, SparseMatrix &RegMat, const std::vector<eslocal> &fixPoints, const Coordinates &coordinates, size_t subdomain)
{
	SparseMatrix Nt; // CSR matice s DOFY
	Nt.rows = 6;
	Nt.cols = K.cols;
	Nt.nnz  = 9 * fixPoints.size();
	Nt.type = 'G';

	std::vector<eslocal> &ROWS = Nt.CSR_I_row_indices;
	std::vector<eslocal> &COLS = Nt.CSR_J_col_indices;
	std::vector<double>  &VALS = Nt.CSR_V_values;

	ROWS.reserve(Nt.rows + 1);
	COLS.reserve(Nt.nnz);
	VALS.reserve(Nt.nnz);

	ROWS.push_back(1);
	ROWS.push_back(ROWS.back() + fixPoints.size());
	ROWS.push_back(ROWS.back() + fixPoints.size());
	ROWS.push_back(ROWS.back() + fixPoints.size());
	ROWS.push_back(ROWS.back() + 2 * fixPoints.size());
	ROWS.push_back(ROWS.back() + 2 * fixPoints.size());
	ROWS.push_back(ROWS.back() + 2 * fixPoints.size());

	for (size_t c = 0; c < 3; c++) {
		std::vector<double> kernel = { 0, 0, 0 };
		kernel[c] = 1;

		for (size_t i = 0; i < fixPoints.size(); i++) {
			COLS.push_back(3 * fixPoints[i] + 1);
			COLS.push_back(3 * fixPoints[i] + 2);
			COLS.push_back(3 * fixPoints[i] + 3);
			VALS.insert(VALS.end(), kernel.begin(), kernel.end());
		}
	}

	for (size_t i = 0; i < fixPoints.size(); i++) {
		const Point &p = coordinates.get(fixPoints[i], subdomain);
		COLS.push_back(3 * fixPoints[i] + 1);
		COLS.push_back(3 * fixPoints[i] + 2);
		VALS.push_back(-p.y);
		VALS.push_back( p.x);
	}

	for (size_t i = 0; i < fixPoints.size(); i++) {
		const Point &p = coordinates.get(fixPoints[i], subdomain);
		COLS.push_back(3 * fixPoints[i] + 1);
		COLS.push_back(3 * fixPoints[i] + 3);
		VALS.push_back(-p.z);
		VALS.push_back( p.x);
	}

	for (size_t i = 0; i < fixPoints.size(); i++) {
		const Point &p = coordinates.get(fixPoints[i], subdomain);
		COLS.push_back(3 * fixPoints[i] + 2);
		COLS.push_back(3 * fixPoints[i] + 3);
		VALS.push_back(-p.z);
		VALS.push_back( p.y);
	}

	SparseMatrix N;
	Nt.MatTranspose( N );

	RegMat.MatMat(Nt, 'N', N);
	RegMat.MatTranspose();
	RegMat.RemoveLower();

	SparseSolverCPU NtN;
	NtN.ImportMatrix(RegMat);
	RegMat.Clear();

	NtN.Factorization("Create RegMat");
	NtN.SolveMat_Sparse(Nt);
	NtN.Clear();

	RegMat.MatMat(N, 'N', Nt);
	RegMat.MatScale(K.getDiagonalMaximum());
}

static void algebraicKernelsAndRegularization(SparseMatrix &K, SparseMatrix &RegMat, SparseMatrix &R, size_t subdomain)
{
	double norm;
	eslocal defect;

	K.get_kernel_from_K(K, RegMat, R, norm, defect, subdomain);
}

void LinearElasticity::composeSubdomain(size_t subdomain)
{
	eslocal subdomainSize = DOFs * _mesh.coordinates().localSize(subdomain);

	SparseVVPMatrix<eslocal> _K;
	DenseMatrix Ke;
	std::vector<double> fe;

	_K.resize(subdomainSize, subdomainSize);
	f[subdomain].resize(subdomainSize);

	const std::vector<eslocal> &partition = _mesh.getPartition();
	const std::vector<Element*> &elements = _mesh.getElements();

	for (eslocal e = partition[subdomain]; e < partition[subdomain + 1]; e++) {

		processElement(Ke, fe, _mesh, subdomain, elements[e]);

		for (size_t i = 0; i < DOFs * elements[e]->size(); i++) {
			size_t row = DOFs * (elements[e]->node(i % elements[e]->size())) + i / elements[e]->size();
			for (size_t j = 0; j < DOFs * elements[e]->size(); j++) {
				size_t column = DOFs * (elements[e]->node(j % elements[e]->size())) + j / elements[e]->size();
				_K(row, column) = Ke(i, j);
			}
			f[subdomain][row] += fe[i];
		}
	}

	// TODO: make it direct
	SparseCSRMatrix<eslocal> csrK = _K;
	K[subdomain] = csrK;

	switch (config::solver::REGULARIZATION) {
	case config::solver::REGULARIZATIONalternative::FIX_POINTS:
		analyticsKernels(R1[subdomain], _mesh.coordinates(), subdomain);
		analyticsRegMat(K[subdomain], RegMat[subdomain], _mesh.computeFixPoints(subdomain, config::mesh::FIX_POINTS), _mesh.coordinates(), subdomain);
		K[subdomain].RemoveLower();
		RegMat[subdomain].RemoveLower();
		K[subdomain].MatAddInPlace(RegMat[subdomain], 'N', 1);
		RegMat[subdomain].ConvertToCOO(1);
		break;
	case config::solver::REGULARIZATIONalternative::NULL_PIVOTS:
		K[subdomain].RemoveLower();
		algebraicKernelsAndRegularization(K[subdomain], RegMat[subdomain], R1[subdomain], subdomain);
		break;
	}

	// TODO:
	R1H[subdomain] = R1[subdomain];
}




