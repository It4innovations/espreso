
#include "assembler.h"

using namespace espreso;

static void processElement(DenseMatrix &Ah, DenseMatrix &B1h, DenseMatrix &B2h, DenseMatrix &Eh, std::vector<double> &fe, const espreso::Mesh &mesh, size_t subdomain, const Element* element)
{
	DenseMatrix coordinates;

	coordinates.resize(element->size(), 2);
	for (size_t i = 0; i < element->size(); i++) {
		coordinates(i, 0) = mesh.coordinates().get(element->node(i), subdomain).x;
		coordinates(i, 1) = mesh.coordinates().get(element->node(i), subdomain).y;
	}

	eslocal Ksize = 3 * element->size();

	fe.resize(Ksize);
	fill(fe.begin(), fe.end(), 0);

	auto node = [ & ] (eslocal n) {
		return mesh.coordinates().get(element->node(n), subdomain);
	};

	double
	x10 = node(1).x - node(0).x,
	x21 = node(2).x - node(1).x,
	x02 = node(0).x - node(2).x,
	y01 = node(0).y - node(1).y,
	y12 = node(1).y - node(2).y,
	y20 = node(2).y - node(0).y,

	tarea = (x10 * y20 - x02 * y01) / 2;

	double alpha = 0, nu = 0.01;
	DenseMatrix xt(3, 1), yt(3, 1), zt(3, 1);
	xt(0, 0) = x21;
	xt(1, 0) = x02;
	xt(2, 0) = x10;

	yt(0, 0) = y12;
	yt(1, 0) = y20;
	yt(2, 0) = y01;

	zt(0, 0) = zt(1, 0) = zt(2, 0) = (3 / 20.0) * alpha * tarea;

	double omega = (81 / 40.0) * nu *(y12 * y12 + x21 * x21 + y20 * y20 + x02 * x02 + y12 * y20 + x21 * x02) / tarea + (81 / 280.0) * alpha * tarea;

	DenseMatrix ytt(3, 3), xtt(3, 3), ztt(3, 3), ytztt, xtztt;
	xtt.multiply(xt, xt, 1, 0, false, true);
	ytt.multiply(yt, yt, 1, 0, false, true);
	ztt.multiply(zt, zt, 1, 0, false, true);
	ytztt.multiply(yt, zt, 1, 0, false, true);
	xtztt.multiply(xt, zt, 1, 0, false, true);

	for (size_t i = 0; i < 3; i++) {
		for (size_t j = 0; j < 3; j ++) {
			Ah(i, j)  = (nu * 0.25 / tarea) * (xtt(i, j) + ytt(i, j)) + ztt(i, j) / omega;
			B1h(i, j) = - (1 / 6.0) * yt(j, 0) - (9.0 / (40 * omega)) * ytztt(i, j);
			B2h(i, j) = - (1 / 6.0) * xt(j, 0) - (9.0 / (40 * omega)) * xtztt(i, j);
			Eh(i, j)  = 81.0 / (1600 * omega) * (ytt(i, j) + xtt(i, j));
		}
	}
}

static void analyticsKernels(SparseMatrix &R1, size_t nodes)
{
	R1.rows = 3 * nodes;
	R1.cols = 2;
	R1.nnz = R1.rows * R1.cols;
	R1.type = 'G';

	R1.dense_values.reserve(R1.nnz);

	std::vector<double> kernel;

	kernel = { 1, 0, 0 };
	for (size_t i = 0; i < nodes; i++) {
		R1.dense_values.insert(R1.dense_values.end(), kernel.begin(), kernel.end());
	}

	kernel = { 0, 1, 0 };
	for (size_t i = 0; i < nodes; i++) {
		R1.dense_values.insert(R1.dense_values.end(), kernel.begin(), kernel.end());
	}
}

static void analyticsRegMat(SparseMatrix &K, SparseMatrix &RegMat)
{
	RegMat.rows = K.rows;
	RegMat.cols = K.cols;
	RegMat.nnz  = 2;
	RegMat.type = K.type;

	RegMat.I_row_indices.push_back(1);
	RegMat.J_col_indices.push_back(1);
	RegMat.I_row_indices.push_back(2);
	RegMat.J_col_indices.push_back(2);
	RegMat.V_values.resize(2, K.getDiagonalMaximum());
	RegMat.ConvertToCSR(1);
}

static void algebraicKernelsAndRegularization(SparseMatrix &K, SparseMatrix &RegMat, SparseMatrix &R, size_t subdomain)
{
	double norm;
	eslocal defect;

	K.get_kernel_from_K(K, RegMat, R, norm, defect, subdomain);
}

void Stokes::composeSubdomain(size_t subdomain)
{
	eslocal subdomainSize = 3 * _mesh.coordinates().localSize(subdomain);
	SparseVVPMatrix<eslocal> _K;
	DenseMatrix Ah(3, 3), B1h(3, 3), B2h(3, 3), Eh(3, 3);
	std::vector<double> fe;

	_K.resize(subdomainSize, subdomainSize);
	f[subdomain].resize(subdomainSize);

	const std::vector<eslocal> &partition = _mesh.getPartition();
	const std::vector<Element*> &elements = _mesh.getElements();

	for (eslocal i = partition[subdomain]; i < partition[subdomain + 1]; i++) {

		processElement(Ah, B1h, B2h, Eh, fe, _mesh, subdomain, elements[i]);

		for (size_t i = 0; i < elements[i]->size(); i++) {
			eslocal row = 3 * elements[i]->node(i);
			for (size_t j = 0; j < elements[i]->size(); j++) {
				eslocal column = 3 * elements[i]->node(j);
				_K(row + 0, column + 0) =   Ah(i, j);
				_K(row + 1, column + 1) =   Ah(i, j);

				_K(row + 2, column + 2) = - Eh(i, j);

				_K(row + 0, column + 2) =   B1h(j, i);
				_K(row + 2, column + 0) =   B1h(i, j);

				_K(row + 1, column + 2) =   B2h(j, i);
				_K(row + 2, column + 1) =   B2h(i, j);
			}
			f[subdomain][row] += fe[i];
		}
	}

	// TODO: make it direct
	SparseCSRMatrix<eslocal> csrK = _K;
	K[subdomain] = csrK;

	switch (config::solver::REGULARIZATION) {
	case config::solver::REGULARIZATIONalternative::FIX_POINTS:
		analyticsKernels(R1[subdomain], _mesh.coordinates().localSize(subdomain));
		analyticsRegMat(K[subdomain], RegMat[subdomain]);
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

	R1H[subdomain] = R1[subdomain];
}




