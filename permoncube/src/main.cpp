#include "permoncube.h"

permoncube::Settings settings;

void setParams(int argc, char** argv)
{
	if (argc != 10) {
		return;
	}

	int _cluster;
	int _subdomains;
	int _elementsInSub;

	for (int i = 0; i < 3; i++) {
		sscanf(argv[i + 1], "%i", &_cluster);
		sscanf(argv[i + 4], "%i", &_subdomains);
		sscanf(argv[i + 7], "%i", &_elementsInSub);
		settings.clusters[i] = _cluster;
		settings.subdomainsInCluster[i] = _subdomains;
		settings.elementsInSubdomain[i] = _elementsInSub;
	}
}

void test_hexa8();
void test_tetra4();
void test_tetra10();

int main(int argc, char** argv)
{
	setParams(argc, argv);
	std::cout << settings;

	test_hexa8();
	test_tetra4();
	test_tetra10();
}

void test_hexa8()
{
	mesh::Mesh mesh;

	size_t cluster[3] = { 0, 0, 0 };
	if (settings.clusters[1] > 1) {
		cluster[1] = 1;
	}
	if (settings.clusters[2] > 2) {
		cluster[2] = 2;
	}
	permoncube::ElementGenerator<permoncube::Hexahedron8> generator(settings);

	generator.mesh(mesh, cluster);

	mesh.saveVTK("hexa8full.vtk", 0.9);
	std::cout << "hexa8full saved\n";

	mesh::SurfaceMesh sMesh(mesh);

	sMesh.saveVTK("hexa8surface.vtk", 0.9);
	std::cout << "hexa8surface saved\n";
}

void test_tetra4()
{
	mesh::Mesh mesh;

	size_t cluster[3] = { 0, 0, 0 };
	if (settings.clusters[1] > 1) {
		cluster[1] = 1;
	}
	if (settings.clusters[2] > 2) {
		cluster[2] = 2;
	}
	permoncube::ElementGenerator<permoncube::Tetrahedron4> generator(settings);

	generator.mesh(mesh, cluster);

	mesh.saveVTK("tetra4full.vtk", 0.9);
	std::cout << "tetra4full saved\n";

	mesh::SurfaceMesh sMesh(mesh);

	sMesh.saveVTK("tetra4surface.vtk", 0.9);
	std::cout << "tetra4surface saved\n";
}

void test_tetra10()
{
	mesh::Mesh mesh;

	size_t cluster[3] = { 0, 0, 0 };
	if (settings.clusters[1] > 1) {
		cluster[1] = 1;
	}
	if (settings.clusters[2] > 2) {
		cluster[2] = 2;
	}
	permoncube::ElementGenerator<permoncube::Tetrahedron10> generator(settings);

	generator.mesh(mesh, cluster);
	std::cout << "generated\n";

	mesh.saveVTK("tetra10full.vtk", 0.9);
	std::cout << "tetra10full saved\n";

	mesh::SurfaceMesh sMesh(mesh);

	sMesh.saveVTK("tetra10surface.vtk", 0.9);
	std::cout << "tetra10surface saved\n";
}
