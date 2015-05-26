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

void test_tetra10();

int main(int argc, char** argv)
{
	setParams(argc, argv);

	permoncube::Settings settings;
	permoncube::Generator *generator = new permoncube::ElementGenerator<permoncube::Tetrahedron4>(settings);

	mesh::Mesh mesh;
	size_t cluster[3] = { 0, 0, 0 };

	std::cout << settings;

	generator->mesh(mesh, cluster);

	//std::cout << mesh.coordinates();

	mesh.saveVTK("permon.vtk", 0.9);

	delete generator;
}


void test_tetra10()
{
	mesh::Mesh m;

	//Permoncube::PM::tetrahedrons10(mesh, mesh.coordinates(), subdomains, elementsInSub);

	int dimension = m.getPartNodesCount(0) * mesh::Point::size();

	SparseCSRMatrix K(dimension, dimension);
	SparseCSRMatrix M(dimension, dimension);
	std::vector<double> f(dimension);

	m.elasticity(K, M, f, 0);

	m.saveVTK("mesh.vtk");

	std::ofstream fileK("K15.txt");
	fileK << K;
	fileK.close();
}
