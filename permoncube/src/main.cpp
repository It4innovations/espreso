#include "permoncube.h"

int subdomains[3] = { 1, 1, 1 };
int elementsInSub[3] = { 5, 5, 5 };

void setParams(int argc, char** argv)
{
	if (argc != 7) {
		return;
	}

	int _subdomains;
	int _elementsInSub;

	for (int i = 0; i < 3; i++) {
		sscanf(argv[i + 1], "%i", &_subdomains);
		sscanf(argv[i + 4], "%i", &_elementsInSub);
		subdomains[i] = _subdomains;
		elementsInSub[i] = _elementsInSub;
	}
}

int main(int argc, char** argv)
{
	setParams(argc, argv);

	Mesh mesh;

	Permoncube::tetrahedrons10(mesh, mesh.coordinates(), subdomains, elementsInSub);

	int dimension = mesh.getPartNodesCount(0) * Point::size();

	SparseCSRMatrix K(dimension, dimension);
	SparseCSRMatrix M(dimension, dimension);
	std::vector<double> f(dimension);

	mesh.elasticity(K, M, f, 0);

	mesh.saveVTK("mesh.vtk");

	std::ofstream fileK("K15.txt");
	fileK << K;
	fileK.close();
}



