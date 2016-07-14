
#include "assembler.h"
//#include "esbem.h"

using namespace espreso;

static void computeKernels(SparseMatrix &R1, const Coordinates &coordinates, size_t subdomain)
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

static void computeRegMat(SparseMatrix &K, SparseMatrix &RegMat, const std::vector<eslocal> &fixPoints, const Coordinates &coordinates, size_t subdomain)
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


void LameSteklovPoincare::composeSubdomain(size_t subdomain)
{
	const std::vector<eslocal> &partition = _mesh.getPartition();
	const std::vector<Element*> &elements = _mesh.getElements();

	DenseMatrix _K;
	eslocal nK = _mesh.coordinates().localSize(subdomain) * Point::size();
	eslocal eSize = partition[subdomain + 1] - partition[subdomain];
	_K.resize(nK, nK);
	std::vector<double> nodes(nK);
	std::vector<eslocal> elems(3 * eSize);

	for (size_t i = 0; i < _mesh.coordinates().localSize(subdomain); i++) {
		nodes[i * Point::size() + 0] = _mesh.coordinates().get(i, subdomain).x;
		nodes[i * Point::size() + 1] = _mesh.coordinates().get(i, subdomain).y;
		nodes[i * Point::size() + 2] = _mesh.coordinates().get(i, subdomain).z;
	}
	for (size_t i = partition[subdomain], index = 0; i < partition[subdomain + 1]; i++, index++) {
		for (size_t j = 0; j < elements[i]->size(); j++) {
			elems[3 * index + j] = elements[i]->node(j);
		}
	}
	ESINFO(GLOBAL_ERROR) << "missing BEM library";
/*
	bem4i::getLameSteklovPoincare<eslocal, double>(
			_K.values(),
			_mesh.coordinates().localSize(subdomain),
			&nodes[0],
			eSize,
			&elems[0],
			0.3,			// nu
			2.1e5,			// E
			3,				// order near
			4,				// order far
			false			// verbose
			);
*/
	DenseMatrix tmp = _K;
	eslocal n = _K.rows();
	for (eslocal i = 0; i < n / 3; i++) {
		for (eslocal j = 0; j < n; j++) {
			tmp(3 * i + 0, j) = _K(0 * (n / 3) + i, j);
			tmp(3 * i + 1, j) = _K(1 * (n / 3) + i, j);
			tmp(3 * i + 2, j) = _K(2 * (n / 3) + i, j);
		}
	}

	for (eslocal i = 0; i < n / 3; i++) {
		for (eslocal j = 0; j < n; j++) {
			_K(j, 3 * i + 0) = tmp(j, 0 * (n / 3) + i);
			_K(j, 3 * i + 1) = tmp(j, 1 * (n / 3) + i);
			_K(j, 3 * i + 2) = tmp(j, 2 * (n / 3) + i);
		}
	}

	// TODO: make it direct
	SparseCSRMatrix<eslocal> csrK = _K;
	K[subdomain] = csrK;

	f[subdomain].clear();
	f[subdomain].resize(_K.rows(), 0);

	computeKernels(R1[subdomain], _mesh.coordinates(), subdomain);
	computeRegMat(K[subdomain], RegMat[subdomain], _mesh.getFixPoints()[subdomain], _mesh.coordinates(), subdomain);

	K[subdomain].RemoveLower();
	RegMat[subdomain].RemoveLower();
}
