#include "esmesh.h"

void test_matrices();
void test_BEM();
void test_meshes();
void test_ansys();
void test_saveData();
void test_ijv_sort();

int main(int argc, char** argv)
{
	test_ijv_sort();
	return 0;
}

void test_ijv_sort()
{
	size_t rows = 50, columns = 50;

	DenseMatrix d(rows, columns);

	for (size_t i = 0; i < rows; i++) {
		for (size_t j = 0; j < columns; j++) {
			if ((j * columns + i) % 15 == 0) {
				d(i, j) = j * columns + i;
			}
		}
	}
	SparseIJVMatrix<eslocal> m(d);
	m.transpose();

	std::cout << m;
}

void test_saveData()
{
	eslocal partsCount = 4;
	eslocal fixPointsCount = 4;

	mesh::Mesh m("matrices/TET/5/elem", "matrices/TET/5/coord", partsCount, fixPointsCount);

	m.saveData();

	//mesh::Mesh m2;
	m.loadData("mesh_0.dat");
	m.saveVTK("part0.vtk");
}

void test_ansys()
{
	eslocal partsCount = 4;
	eslocal fixPointsCount = 4;

	mesh::Ansys ansys("matrices/spanner/Model");
	ansys.coordinatesProperty(mesh::CP::DIRICHLET_X) = "BC/Elasticity/NUX.dat";
	ansys.coordinatesProperty(mesh::CP::DIRICHLET_Y) = "BC/Elasticity/NUY.dat";
	ansys.coordinatesProperty(mesh::CP::DIRICHLET_Z) = "BC/Elasticity/NUZ.dat";

	std::cout << ansys;

	mesh::Mesh m(ansys, partsCount, fixPointsCount);

	mesh::Boundaries b(m);

	mesh::SurfaceMesh sMesh(m);

	mesh::CommonFacesMesh cMesh(sMesh);

	mesh::CornerLinesMesh lMesh(m);

	m.saveVTK("mesh.vtk", 0.6);
	sMesh.saveVTK("surface.vtk", 0.6);
	cMesh.saveVTK("faces.vtk", 0.6);
	lMesh.saveVTK("lines.vtk", 0.6);
}


void test_meshes()
{
	eslocal partsCount = 4;
	eslocal fixPointsCount = 4;

	mesh::Mesh m("matrices/TET/10/elem", "matrices/TET/10/coord", partsCount, fixPointsCount);

	mesh::Boundaries b(m);

	mesh::SurfaceMesh sMesh(m);

	mesh::CommonFacesMesh cMesh(sMesh);

	mesh::CornerLinesMesh lMesh(m);

	m.saveVTK("mesh.vtk", 0.6);
	sMesh.saveVTK("surface.vtk", 0.6);
	cMesh.saveVTK("faces.vtk", 0.6);
	lMesh.saveVTK("lines.vtk", 0.6);
}

void test_BEM()
{
	eslocal partsCount = 4;
	eslocal fixPointsCount = 4;

	mesh::Mesh m("matrices/TET/10/elem", "matrices/TET/10/coord", partsCount, fixPointsCount);

	mesh::SurfaceMesh bMesh;
	m.getSurface(bMesh);

	std::vector<DenseMatrix> K_mat;

	K_mat.reserve(partsCount);
	for (eslocal d = 0; d < partsCount; d++) {
		K_mat.push_back( DenseMatrix (0, 0) );
	}

	for (eslocal d = 0; d < partsCount; d++) {

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
	SparseDOKMatrix<eslocal> dok;
	SparseVVPMatrix<eslocal> vvp;

	matrices.push_back(&d);
	fill_matrix(matrices.back());

	matrices.push_back(&dok);
	fill_matrix(matrices.back());

	matrices.push_back(&vvp);
	fill_matrix(matrices.back());
	matrices.back()->set(1, 1, 3);
	matrices.back()->set(1, 1, -3);

	// IJV matrix
	matrices.push_back(new SparseIJVMatrix<eslocal>(d));
	matrices.push_back(new SparseIJVMatrix<eslocal>());
	*dynamic_cast<SparseIJVMatrix<eslocal>*>(matrices.back()) = d;
	matrices.push_back(new SparseIJVMatrix<eslocal>(dok));
	matrices.push_back(new SparseIJVMatrix<eslocal>());
	*dynamic_cast<SparseIJVMatrix<eslocal>*>(matrices.back()) = dok;
	matrices.push_back(new SparseIJVMatrix<eslocal>(vvp));
	matrices.push_back(new SparseIJVMatrix<eslocal>());
	*dynamic_cast<SparseIJVMatrix<eslocal>*>(matrices.back()) = vvp;

	// CSR matrix
	matrices.push_back(new SparseCSRMatrix<eslocal>(d));
	matrices.push_back(new SparseCSRMatrix<eslocal>());
	*dynamic_cast<SparseCSRMatrix<eslocal>*>(matrices.back()) = d;
	matrices.push_back(new SparseCSRMatrix<eslocal>(dok));
	matrices.push_back(new SparseCSRMatrix<eslocal>());
	*dynamic_cast<SparseCSRMatrix<eslocal>*>(matrices.back()) = dok;

	SparseCSRMatrix<eslocal> ccc(d);
	SparseIJVMatrix<eslocal> iii(d);

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

	SparseCSRMatrix<eslocal> csr(d);
	SparseCSRMatrix<eslocal> csrT(csr);
	csrT.transpose();
	SparseCSRMatrix<eslocal> csrTT(csrT);
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

	SparseIJVMatrix<eslocal> ijv(d);
	SparseIJVMatrix<eslocal> ijvT(csr);
	ijvT.transpose();

	for (size_t r = 0; r < ijvT.rows(); r++) {
		for (size_t c = 0; c < ijvT.columns(); c++) {
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

	SparseDOKMatrix<eslocal> dokA(m, k);
	SparseDOKMatrix<eslocal> dokB(k, n);
	SparseDOKMatrix<eslocal> dokResult(m, n);
	eslocal result[] = { 1, 4, 10, 20, 35 };

	for (eslocal i = 0; i < m; i++) {
		for (eslocal j = i; j < n; j++) {
			dokResult(i, j) = result[j - i];
		}
	}

	for (eslocal i = 0; i < m; i++) {
		for (eslocal j = 0; j < k - i; j++) {
			dokA(i, j + i) = j + 1;
		}
	}

	for (eslocal i = 0; i < k; i++) {
		for (eslocal j = 0; j < n - i; j++) {
			dokB(i, j + i) = j + 1;
		}
	}

	SparseCSRMatrix<eslocal> A(dokA);
	SparseCSRMatrix<eslocal> B(dokB);
	SparseCSRMatrix<eslocal> C;

	C.multiply(A, B);
	DenseMatrix denseC(C);

	for (eslocal i = 0; i < m; i++) {
		for (eslocal j = 0; j < n; j++) {
			if (denseC(i, j) != dokResult(i, j)) {
				std::cerr << "CSR A * CSR B is incorrect\n";
				exit(EXIT_FAILURE);
			}
		}
	}

	A.transpose();

	SparseCSRMatrix<eslocal> D;

	D.multiply(A, B, true);

	for (eslocal i = 0; i < m; i++) {
		for (eslocal j = 0; j < n; j++) {
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

	for (eslocal i = 0; i < m; i++) {
		for (eslocal j = 0; j < n; j++) {
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

	for (eslocal i = 0; i < m; i++) {
		for (eslocal j = 0; j < n; j++) {
			if (dAB.get(i, j) != dokResult(i, j)) {
				std::cerr << "dense: AT * B is incorrect\n";
				exit(EXIT_FAILURE);
			}
		}
	}

	dAB.multiply(dA, dBT, 1, 0, false, true);

	for (eslocal i = 0; i < m; i++) {
		for (eslocal j = 0; j < n; j++) {
			if (dAB.get(i, j) != dokResult(i, j)) {
				std::cerr << "dense: A * BT is incorrect\n";
				exit(EXIT_FAILURE);
			}
		}
	}

	dAB.multiply(dAT, dBT, 1, 0, true, true);

	for (eslocal i = 0; i < m; i++) {
		for (eslocal j = 0; j < n; j++) {
			if (dAB.get(i, j) != dokResult(i, j)) {
				std::cerr << "dense: AT * BT is incorrect\n";
				exit(EXIT_FAILURE);
			}
		}
	}

	dAB.multiply(dAT, dBT, 1, 1, true, true);

	for (eslocal i = 0; i < m; i++) {
		for (eslocal j = 0; j < n; j++) {
			if (dAB.get(i, j) != 2 * dokResult(i, j)) {
				std::cerr << "dense: AT * BT + C  is incorrect\n";
				exit(EXIT_FAILURE);
			}
		}
	}

}



