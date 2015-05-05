#include "esmesh.h"

int main(int argc, char** argv)
{
	int partsCount = 4;
	int fixPointsCount = 8;

	Coordinates coords("matrices/HEX/15/coord");
	Mesh mesh(coords);
	mesh = Mesh("matrices/HEX/15/elem", coords, partsCount, fixPointsCount);

	int dimension = mesh.getPartNodesCount(0) * Point::size();

	SparseCSRMatrix K(dimension, dimension);
	SparseCSRMatrix M(dimension, dimension);
	std::vector<double> f(dimension);

	mesh.assemble_matrix(K, M, f, 0);

	//mesh.saveVTK();
}
