#include "esmesh.h"

int main(int argc, char** argv)
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
