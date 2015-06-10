
#include <set>

#include "espmcube.h"


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

void test_hexa8(size_t cluster[]);
void test_hexa20(size_t cluster[]);
void test_tetra4(size_t cluster[]);
void test_tetra10(size_t cluster[]);
void test_boudaries();

int main(int argc, char** argv)
{
	setParams(argc, argv);

	size_t cluster[3] = { 0, 0, 0 };
	if (settings.clusters[1] > 1) {
		cluster[1] = 1;
	}
	if (settings.clusters[2] > 2) {
		cluster[2] = 2;
	}

	//test_boudaries();

	std::cout << settings;
	//test_hexa8(cluster);
	test_hexa20(cluster);
	//test_tetra4(cluster);
	//test_tetra10(cluster);
}

void test_boudaries()
{
	for (int i = 0; i < 3; i++) {
		settings.clusters[i] = 3;
		settings.subdomainsInCluster[i] = 5;
		settings.elementsInSubdomain[i] = 1; // Mandatory !!!
	}
	std::cout << settings;

	// MESH BOUNDARIES
	mesh::Mesh mesh;
	permoncube::Settings _settings = settings;
	for (int i = 0; i < 3; i++) {
		_settings.clusters[i] = 1;
		_settings.subdomainsInCluster[i] = settings.clusters[i];
		_settings.elementsInSubdomain[i] = settings.subdomainsInCluster[i];
	}
	permoncube::ElementGenerator<permoncube::Tetrahedron10> g1(_settings);

	size_t cluster[3] = { 0, 0, 0 };
	g1.mesh(mesh, cluster);
	mesh::Boundaries bMesh(mesh);

	// GLOBAL BOUNDARIES
	permoncube::ElementGenerator<permoncube::Tetrahedron10> g2(settings);
	mesh::Boundaries bGlobal;
	g2.fillGlobalBoundaries(bGlobal);

	if (bMesh.size() != bGlobal.size()) {
		std::cerr << "Incorrect boundaries size in Permoncube\n";
		std::cerr << bMesh.size() << " vs. " << bGlobal.size() << "\n";
		exit(EXIT_FAILURE);
	}
	for (size_t i = 0; i < bMesh.size(); i++) {
		if (bMesh[i].size() != bGlobal[i].size()) {
			std::cerr << "Incorrect boundaries node size in Permoncube\n";
			exit(EXIT_FAILURE);
		}
		std::set<eslocal>::const_iterator it;
		for (it = bMesh[i].begin(); it != bMesh[i].end(); ++it) {
			if (!bGlobal[i].count(*it)) {
				std::cerr << "Incorrect boundaries node occurrence in Permoncube\n";
				exit(EXIT_FAILURE);
			}
		}
	}
}

void test_hexa8(size_t cluster[])
{
	mesh::Mesh mesh;

	permoncube::ElementGenerator<permoncube::Hexahedron8> generator(settings);

	generator.mesh(mesh, cluster);

	mesh.saveVTK("hexa8full.vtk", 0.9);
	std::cout << "hexa8full saved\n";

	mesh::SurfaceMesh sMesh(mesh);

	sMesh.saveVTK("hexa8surface.vtk", 0.9);
	std::cout << "hexa8surface saved\n";
}

void test_hexa20(size_t cluster[])
{
	mesh::Mesh mesh;

	permoncube::ElementGenerator<permoncube::Hexahedron20> generator(settings);

	generator.mesh(mesh, cluster);

	mesh.saveVTK("hexa20full.vtk", 0.9);
	std::cout << "hexa20full saved\n";

	mesh::SurfaceMesh sMesh(mesh);

	sMesh.saveVTK("hexa20surface.vtk", 0.9);
	std::cout << "hexa20surface saved\n";
}

void test_tetra4(size_t cluster[])
{
	mesh::Mesh mesh;

	permoncube::ElementGenerator<permoncube::Tetrahedron4> generator(settings);

	generator.mesh(mesh, cluster);

	mesh.saveVTK("tetra4full.vtk", 0.9);
	std::cout << "tetra4full saved\n";

	mesh::SurfaceMesh sMesh(mesh);

	sMesh.saveVTK("tetra4surface.vtk", 0.9);
	std::cout << "tetra4surface saved\n";
}

void test_tetra10(size_t cluster[])
{
	mesh::Mesh mesh;

	permoncube::ElementGenerator<permoncube::Tetrahedron10> generator(settings);

	generator.mesh(mesh, cluster);
	std::cout << "generated\n";

	mesh.saveVTK("tetra10full.vtk", 0.9);
	std::cout << "tetra10full saved\n";

	mesh::SurfaceMesh sMesh(mesh);

	sMesh.saveVTK("tetra10surface.vtk", 0.9);
	std::cout << "tetra10surface saved\n";
}
