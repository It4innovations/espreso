#include "esmesh.h"

int main(int argc, char** argv)
{
	int partsCount = 4;
	int fixPointsCount = 8;

	Coordinates coords("matrices/HEX/5/coord");
	Mesh mesh(coords);
	mesh = Mesh("matrices/HEX/5/elem", coords, partsCount, fixPointsCount);

	int dimension = mesh.getPartNodesCount(0) * Point::size();

	SparseCSRMatrix K(dimension, dimension);
	SparseCSRMatrix M(dimension, dimension);
	std::vector<double> f(dimension);

	//mesh.assemble_matrix(K, M, f, 0);

	mesh.saveVTK("mesh.vtk");

	Coordinates c;
	Mesh bem(c);
	mesh.getBEM(bem);

	bem.saveVTK("bem.vtk");

	//mesh.saveVTK();
}
