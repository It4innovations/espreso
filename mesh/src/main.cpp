#include "esmesh.h"

void test_matrices();
void test_BEM();
void test_meshes();

int main(int argc, char** argv)
{
	test_meshes();
	return 0;
}

void test_meshes()
{
	int partsCount = 4;
	int fixPointsCount = 4;

	mesh::Mesh m("matrices/HEX/5/elem", "matrices/HEX/5/coord", partsCount, fixPointsCount);

	mesh::Boundaries b(m);

	mesh::SurfaceMesh sMesh(m);

	mesh::CommonFacesMesh cMesh(sMesh);

	m.saveVTK("mesh.vtk", 0.6);
	sMesh.saveVTK("surface.vtk", 0.6);
	cMesh.saveVTK("faces.vtk", 0.6);
}

void test_BEM()
{
	int partsCount = 4;
	int fixPointsCount = 4;

	mesh::Mesh m("matrices/TET/10/elem", "matrices/TET/10/coord", partsCount, fixPointsCount);

	mesh::SurfaceMesh bMesh;
	m.getSurface(bMesh);

	std::vector<DenseMatrix> K_mat;

	K_mat.reserve(partsCount);
	for (int d = 0; d < partsCount; d++) {
		K_mat.push_back( DenseMatrix (0, 0) );
	}

	for (int d = 0; d < partsCount; d++) {

		bMesh.elasticity(K_mat[d], d);

		std::cout << d << " " << std::endl;
	}

	bMesh.saveVTK("bem.vtk");

}

void fill_matrix(Matrix *m)
{
	m->resize(5, 4);
	m->operator ()(0, 0) = 1;
	m->operator ()(1, 1) = 2;
	m->operator ()(2, 2) = 3;
	m->operator ()(3, 3) = 4;
	m->operator ()(0, 3) = 3.3;
	m->set(4, 3, 5.5);
}

void test_matrices()
{
	std::vector<Matrix*> matrices;

	DenseMatrix d;
	SparseDOKMatrix dok;
	SparseVVPMatrix vvp;

	matrices.push_back(&d);
	fill_matrix(matrices.back());

	matrices.push_back(&dok);
	fill_matrix(matrices.back());

	matrices.push_back(&vvp);
	fill_matrix(matrices.back());
	matrices.back()->set(1, 1, 3);
	matrices.back()->set(1, 1, -3);

	// IJV matrix
	matrices.push_back(new SparseIJVMatrix(d));
	matrices.push_back(new SparseIJVMatrix());
	*dynamic_cast<SparseIJVMatrix*>(matrices.back()) = d;
	matrices.push_back(new SparseIJVMatrix(dok));
	matrices.push_back(new SparseIJVMatrix());
	*dynamic_cast<SparseIJVMatrix*>(matrices.back()) = dok;
	matrices.push_back(new SparseIJVMatrix(vvp));
	matrices.push_back(new SparseIJVMatrix());
	*dynamic_cast<SparseIJVMatrix*>(matrices.back()) = vvp;

	// CSR matrix
	matrices.push_back(new SparseCSRMatrix(d));
	matrices.push_back(new SparseCSRMatrix());
	*dynamic_cast<SparseCSRMatrix*>(matrices.back()) = d;
	matrices.push_back(new SparseCSRMatrix(dok));
	matrices.push_back(new SparseCSRMatrix());
	*dynamic_cast<SparseCSRMatrix*>(matrices.back()) = dok;

	SparseCSRMatrix ccc(d);
	SparseIJVMatrix iii(d);

	matrices.push_back(new DenseMatrix(ccc));
	matrices.push_back(new DenseMatrix(iii));

	vvp.shrink();
	for (size_t i = 1; i < matrices.size(); i++) {
		for (size_t r = 0; r < matrices[i]->rows(); r++) {
			for (size_t c = 0; c < matrices[i]->columns(); c++) {
				const Matrix *m = matrices[i];
				if (matrices[0]->get(r, c) != m->operator ()(r, c)) {
					std::cerr << *matrices[0];
					std::cerr << *matrices[i];
					std::cerr << "Matrix: " << i << ", ";
					std::cerr << "row: " << r << ", column: " << c << " -> ";
					std::cerr << matrices[0]->get(r, c) << " != ";
					std::cerr << m->operator ()(r, c) << "\n";
					return;
				}
			}
		}
	}

	for (size_t i = 3; i < matrices.size(); i++) {
		delete matrices[i];
	}

	DenseMatrix dT(d);
	dT.transpose();

	for (size_t r = 0; r < d.rows(); r++) {
		for (size_t c = 0; c < d.columns(); c++) {
			if (d(r, c) != dT(c, r)) {
				std::cerr << d;
				std::cerr << dT;
				std::cerr << "Transpose of dense matrix is incorrect.\n";
				return;
			}
		}
	}

	DenseMatrix dokT(dok);
	dokT.transpose();

	for (size_t r = 0; r < d.rows(); r++) {
		for (size_t c = 0; c < d.columns(); c++) {
			if (dok(r, c) != dokT(c, r)) {
				std::cerr << dok;
				std::cerr << dokT;
				std::cerr << "Transpose of DOK matrix is incorrect.\n";
				return;
			}
		}
	}

	SparseCSRMatrix csr(d);
	SparseCSRMatrix csrT(csr);
	csrT.transpose();
	SparseCSRMatrix csrTT(csrT);
	csrTT.transpose();

	for (size_t r = 0; r < d.rows(); r++) {
		for (size_t c = 0; c < d.columns(); c++) {
			if (csr.get(r, c) != csrT.get(c, r)) {
				std::cerr << csr;
				std::cerr << csrT;
				std::cerr << "Transpose of CSR matrix is incorrect.\n";
				return;
			}
		}
	}

	for (size_t r = 0; r < d.rows(); r++) {
		for (size_t c = 0; c < d.columns(); c++) {
			if (csr.get(r, c) != csrTT.get(r, c)) {
				std::cerr << csr;
				std::cerr << csrTT;
				std::cerr << "Transpose of CSR matrix is incorrect.\n";
				return;
			}
		}
	}

	SparseIJVMatrix ijv(d);
	SparseIJVMatrix ijvT(csr);
	ijvT.transpose();

	for (size_t r = 0; r < d.rows(); r++) {
		for (size_t c = 0; c < d.columns(); c++) {
			if (ijv.get(r, c) != ijvT.get(c, r)) {
				std::cerr << ijv;
				std::cerr << ijvT;
				std::cerr << "Transpose of IJV matrix is incorrect.\n";
				return;
			}
		}
	}

	size_t m = 3;
	size_t n = 5;
	size_t k = 6;

	SparseDOKMatrix dokA(m, k);
	SparseDOKMatrix dokB(k, n);
	SparseDOKMatrix dokResult(m, n);
	int result[] = { 1, 4, 10, 20, 35 };

	for (int i = 0; i < m; i++) {
		for (int j = i; j < n; j++) {
			dokResult(i, j) = result[j - i];
		}
	}

	for (int i = 0; i < m; i++) {
		for (int j = 0; j < k - i; j++) {
			dokA(i, j + i) = j + 1;
		}
	}

	for (int i = 0; i < k; i++) {
		for (int j = 0; j < n - i; j++) {
			dokB(i, j + i) = j + 1;
		}
	}

	SparseCSRMatrix A(dokA);
	SparseCSRMatrix B(dokB);
	SparseCSRMatrix C;

	C.multiply(A, B);
	DenseMatrix denseC(C);

	for (int i = 0; i < m; i++) {
		for (int j = 0; j < n; j++) {
			if (denseC(i, j) != dokResult(i, j)) {
				std::cerr << "CSR A * CSR B is incorrect\n";
				exit(EXIT_FAILURE);
			}
		}
	}

	A.transpose();

	SparseCSRMatrix D;

	D.multiply(A, B, true);

	for (int i = 0; i < m; i++) {
		for (int j = 0; j < n; j++) {
			if (D.get(i, j) != dokResult(i, j)) {
				std::cerr << "trans CSR A * CSR B is incorrect\n";
				exit(EXIT_FAILURE);
			}
		}
	}

	DenseMatrix dA(dokA);
	DenseMatrix dAT(dA);
	dAT.transpose();
	DenseMatrix dB(dokB);
	DenseMatrix dBT(dB);
	dBT.transpose();
	DenseMatrix dAB;

	dAB.multiply(dA, dB);

	for (int i = 0; i < m; i++) {
		for (int j = 0; j < n; j++) {
			if (dAB.get(i, j) != dokResult(i, j)) {
				std::cerr << "dense:  A * B is incorrect\n";
				exit(EXIT_FAILURE);
			}
		}
	}

	DenseMatrix dd1(3, 2);
	DenseMatrix dd2(3, 3);

	dd1(0, 0) = 1;
	dd1(1, 0) = 1;
	dd2(1, 0) = 1;
	dd2(1, 1) = 2;
	dd2(1, 2) = 3;

	dAB.multiply(dAT, dB, 1, 0, true);

	for (int i = 0; i < m; i++) {
		for (int j = 0; j < n; j++) {
			if (dAB.get(i, j) != dokResult(i, j)) {
				std::cerr << "dense: AT * B is incorrect\n";
				exit(EXIT_FAILURE);
			}
		}
	}

	dAB.multiply(dA, dBT, 1, 0, false, true);

	for (int i = 0; i < m; i++) {
		for (int j = 0; j < n; j++) {
			if (dAB.get(i, j) != dokResult(i, j)) {
				std::cerr << "dense: A * BT is incorrect\n";
				exit(EXIT_FAILURE);
			}
		}
	}

	dAB.multiply(dAT, dBT, 1, 0, true, true);

	for (int i = 0; i < m; i++) {
		for (int j = 0; j < n; j++) {
			if (dAB.get(i, j) != dokResult(i, j)) {
				std::cerr << "dense: AT * BT is incorrect\n";
				exit(EXIT_FAILURE);
			}
		}
	}

	dAB.multiply(dAT, dBT, 1, 1, true, true);

	for (int i = 0; i < m; i++) {
		for (int j = 0; j < n; j++) {
			if (dAB.get(i, j) != 2 * dokResult(i, j)) {
				std::cerr << "dense: AT * BT + C  is incorrect\n";
				exit(EXIT_FAILURE);
			}
		}
	}

}



