
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
		for (size_t i = 0; i < mesh.fixPoints(p).size(); i++) {
			fixPoints[p].push_back(mesh.fixPoints(p)[i]->node(0));
		}
	}
	output.store(fixPoints, shrinkSubdomain, shringCluster);
}

void VTK_Full::corners(const Mesh &mesh, const std::string &path, double shrinkSubdomain, double shringCluster)
{
	VTK_Full output(mesh, path);
	std::vector<std::vector<eslocal> > corners(mesh.parts());

	for (size_t i = 0; i < mesh.corners().size(); i++) {
		for (size_t d = 0; d < mesh.corners()[i]->domains().size(); d++) {
			corners[mesh.corners()[i]->domains()[d]].push_back(mesh.corners()[i]->node(0));
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

	ESINFO(GLOBAL_ERROR) << "Broken VTK out dirichlet";
//	const std::vector<BoundaryCondition*> &bc = mesh.boundaryConditions();
//	for (size_t i = 0; i < bc.size(); i++) {
//		if (bc[i]->type() == ConditionType::DIRICHLET) {
//			eslocal DOFs = 3;
//			for (size_t p = 0; p < mesh.parts(); p++) {
//				auto &l2c = mesh.coordinates().localToCluster(p);
//				for (size_t j = 0; j < bc[i]->DOFs().size(); j++) {
//					if (bc[i]->DOFs()[j] % DOFs == 0 && std::binary_search(l2c.begin(), l2c.end(), bc[i]->DOFs()[j] / DOFs)) {
//						dx[p].push_back(std::lower_bound(l2c.begin(), l2c.end(), bc[i]->DOFs()[j] / DOFs) - l2c.begin());
//					}
//					if (bc[i]->DOFs()[j] % DOFs == 1 && std::binary_search(l2c.begin(), l2c.end(), bc[i]->DOFs()[j] / DOFs)) {
//						dy[p].push_back(std::lower_bound(l2c.begin(), l2c.end(), bc[i]->DOFs()[j] / DOFs) - l2c.begin());
//					}
//					if (bc[i]->DOFs()[j] % DOFs == 2 && std::binary_search(l2c.begin(), l2c.end(), bc[i]->DOFs()[j] / DOFs)) {
//						dz[p].push_back(std::lower_bound(l2c.begin(), l2c.end(), bc[i]->DOFs()[j] / DOFs) - l2c.begin());
//					}
//				}
//			}
//		}
//	}

	outputx.store(dx, shrinkSubdomain, shringCluster);
	outputy.store(dy, shrinkSubdomain, shringCluster);
	outputz.store(dz, shrinkSubdomain, shringCluster);
}

void VTK_Full::averaging(const Mesh &mesh, const std::string &path, double shrinkSubdomain, double shringCluster)
{
	VTK_Full output(mesh, path);
	std::vector<std::vector<eslocal> > averaging(mesh.parts());

//	for (size_t p = 0; p < mesh.parts(); p++) {
//		for (size_t i = 0; i < mesh.coordinates().localToCluster(p).size(); i++) {
//			if (mesh.subdomainBoundaries().isAveraging(mesh.coordinates().localToCluster(p)[i])) {
//				auto &nodes = mesh.subdomainBoundaries().averaging(mesh.coordinates().localToCluster(p)[i]);
//				for (size_t n = 0; n < nodes.size(); n++) {
//					averaging[p].push_back(mesh.coordinates().localIndex(nodes[n], p));
//				}
//			}
//		}
//	}

	output.store(averaging, shrinkSubdomain, shringCluster);
}



