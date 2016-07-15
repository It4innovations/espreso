
#include "vtk.h"
#include "esconfig.h"

using namespace espreso::output;

void VTK_Full::coordinatesDisplacement(const std::vector<std::vector<double> > &displacement, size_t dofs)
{

	size_t size = 0;
	for (size_t p = 0; p < _mesh.parts(); p++) {
		size += _mesh.coordinates().localSize(p);
	}

	_vtk << "\n";
	_vtk << "POINT_DATA " << size << "\n";
	_vtk << "SCALARS displacements float " << dofs << "\n";
	_vtk << "LOOKUP_TABLE default\n";
	for (size_t p = 0; p < displacement.size(); p++) {
		for (size_t i = 0; i < displacement[p].size() / dofs; i++) {
			for (size_t d = 0; d < dofs; d++) {
				_vtk << displacement[p][dofs * i + d] << " ";
			}
			_vtk << "\n";
		}
	}
}

void VTK_Full::mesh(const Mesh &mesh, const std::string &path, double shrinkSubdomain, double shringCluster)
{
	VTK_Full output(mesh, path);
	output.store(shrinkSubdomain, shringCluster);
}

void VTK_Full::fixPoints(const Mesh &mesh, const std::string &path, double shrinkSubdomain, double shringCluster)
{
	VTK_Full output(mesh, path);
	std::vector<std::vector<eslocal> > fixPoints(mesh.parts());
	for (size_t p = 0; p < mesh.parts(); p++) {
		fixPoints[p] = mesh.computeFixPoints(p, config::mesh::FIX_POINTS);
	}
	output.store(fixPoints, shrinkSubdomain, shringCluster);
}

void VTK_Full::corners(const Mesh &mesh, const std::string &path, double shrinkSubdomain, double shringCluster)
{
	VTK_Full output(mesh, path);
	std::vector<std::vector<eslocal> > corners(mesh.parts());

	for (size_t p = 0; p < mesh.parts(); p++) {
		for (size_t i = 0; i < mesh.coordinates().localToCluster(p).size(); i++) {
			if (mesh.subdomainBoundaries().isCorner(mesh.coordinates().localToCluster(p)[i])) {
				corners[p].push_back(i);
			}
		}
	}

	output.store(corners, shrinkSubdomain, shringCluster);
}

void VTK_Full::dirichlet(const Mesh &mesh, const std::string &path, double shrinkSubdomain, double shringCluster)
{
	VTK_Full outputx(mesh, path + "X");
	VTK_Full outputy(mesh, path + "Y");
	VTK_Full outputz(mesh, path + "Z");

	std::vector<std::vector<eslocal> > dx(mesh.parts());
	std::vector<std::vector<eslocal> > dy(mesh.parts());
	std::vector<std::vector<eslocal> > dz(mesh.parts());

	auto &dxMap = mesh.coordinates().property(DIRICHLET_X).values();
	auto &dyMap = mesh.coordinates().property(DIRICHLET_Y).values();
	auto &dzMap = mesh.coordinates().property(DIRICHLET_Z).values();

	for (size_t p = 0; p < mesh.parts(); p++) {
		auto &l2c = mesh.coordinates().localToCluster(p);
		for (size_t i = 0; i < l2c.size(); i++) {
			if (dxMap.find(l2c[i]) != dxMap.end()) {
				dx[p].push_back(i);
			}
			if (dyMap.find(l2c[i]) != dyMap.end()) {
				dy[p].push_back(i);
			}
			if (dzMap.find(l2c[i]) != dzMap.end()) {
				dz[p].push_back(i);
			}
		}
	}

	const std::vector<BoundaryCondition*> &bc = mesh.boundaryConditions();
	for (size_t i = 0; i < bc.size(); i++) {
		if (bc[i]->type() == ConditionType::DIRICHLET) {
			for (size_t p = 0; p < mesh.parts(); p++) {
				auto &l2c = mesh.coordinates().localToCluster(p);
				for (size_t j = 0; j < bc[i]->DOFs().size(); j++) {
					if (bc[i]->DOFs()[j] % mesh.DOFs() == 0 && std::binary_search(l2c.begin(), l2c.end(), bc[i]->DOFs()[j] / mesh.DOFs())) {
						dx[p].push_back(std::lower_bound(l2c.begin(), l2c.end(), bc[i]->DOFs()[j] / mesh.DOFs()) - l2c.begin());
					}
					if (bc[i]->DOFs()[j] % mesh.DOFs() == 1 && std::binary_search(l2c.begin(), l2c.end(), bc[i]->DOFs()[j] / mesh.DOFs())) {
						dy[p].push_back(std::lower_bound(l2c.begin(), l2c.end(), bc[i]->DOFs()[j] / mesh.DOFs()) - l2c.begin());
					}
					if (bc[i]->DOFs()[j] % mesh.DOFs() == 2 && std::binary_search(l2c.begin(), l2c.end(), bc[i]->DOFs()[j] / mesh.DOFs())) {
						dz[p].push_back(std::lower_bound(l2c.begin(), l2c.end(), bc[i]->DOFs()[j] / mesh.DOFs()) - l2c.begin());
					}
				}
			}
		}
	}

	outputx.store(dx, shrinkSubdomain, shringCluster);
	outputy.store(dy, shrinkSubdomain, shringCluster);
	outputz.store(dz, shrinkSubdomain, shringCluster);
}

void VTK_Full::averaging(const Mesh &mesh, const std::string &path, double shrinkSubdomain, double shringCluster)
{
	VTK_Full output(mesh, path);
	std::vector<std::vector<eslocal> > averaging(mesh.parts());

	for (size_t p = 0; p < mesh.parts(); p++) {
		for (size_t i = 0; i < mesh.coordinates().localToCluster(p).size(); i++) {
			if (mesh.subdomainBoundaries().isAveraging(mesh.coordinates().localToCluster(p)[i])) {
				auto &nodes = mesh.subdomainBoundaries().averaging(mesh.coordinates().localToCluster(p)[i]);
				for (size_t n = 0; n < nodes.size(); n++) {
					averaging[p].push_back(mesh.coordinates().localIndex(nodes[n], p));
				}
			}
		}
	}

	output.store(averaging, shrinkSubdomain, shringCluster);
}



