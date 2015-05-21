#include "esmesh.h"

void test_matrices();
void test_BEM();

int main(int argc, char** argv)
{
	test_matrices();
	return 0;
}

void test_BEM()
{
	int partsCount = 4;
	int fixPointsCount = 4;

	Mesh mesh("matrices/TET/5/elem", "matrices/TET/5/coord", partsCount, fixPointsCount);

	int dimension = mesh.getPartNodesCount(0) * Point::size();

	SparseCSRMatrix K(dimension, dimension);
	SparseCSRMatrix M(dimension, dimension);
	std::vector<double> f(dimension);

	mesh.elasticity(K, M, f, 0);

	mesh.saveVTK("mesh.vtk");

	Boundaries b(mesh);

	BoundaryMesh bMesh;
	mesh.getBoundary(bMesh);

	std::vector<DenseMatrix> K_mat;

	K_mat.reserve(partsCount);
	for (int d = 0; d < partsCount; d++) {
		K_mat.push_back( DenseMatrix (0, 0) );
	}

	for (int d = 0; d < partsCount; d++) {

		bMesh.elasticity(K_mat[d], d);

		std::cout << d << " " << std::endl;
	}

	//bem.saveVTK("bem.vtk");

	//mesh.saveVTK();
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

}



