#include "esmesh.h"

int main(int argc, char** argv)
{
	int partsCount = 1;
	int fixPointsCount = 8;

	Coordinates coords("matrices/TET/10/coord");
	Mesh mesh(coords);
	mesh = Mesh("matrices/TET/10/elem", coords, partsCount, fixPointsCount);

	int dimension = mesh.getPartNodesCount(0) * Point::size();

	SparseCSRMatrix K(dimension, dimension);
	SparseCSRMatrix M(dimension, dimension);
	std::vector<double> f(dimension);

	mesh.elasticity(K, M, f, 0);

	mesh.saveVTK("mesh.vtk");

	Boundaries b(mesh, coords);

	Coordinates c;
	BoundaryMesh bem(c);
	mesh.getBoundary(bem);

	DenseMatrix BK(0, 0);
	bem.elasticity(BK, 0);

	bem.saveVTK("bem.vtk");

	//mesh.saveVTK();
}
