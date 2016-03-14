
#ifndef OUTPUT_VTK_FULL_VTK_H_
#define OUTPUT_VTK_FULL_VTK_H_

#include "../vtk.h"

namespace espreso {
namespace output {

class VTK_Full: public VTK {

public:
	VTK_Full(const Mesh &mesh, const std::string &path): VTK(mesh, path) { };

	static void mesh(const Mesh &mesh, const std::string &path, double shrinkSubdomain, double shringCluster)
	{
		VTK_Full output(mesh, path);
		output.store(shrinkSubdomain, shringCluster);
	}

	static void fixPoints(const Mesh &mesh, const std::string &path, double shrinkSubdomain, double shringCluster)
	{
		VTK_Full output(mesh, path);
		output.store(mesh.getFixPoints(), shrinkSubdomain, shringCluster);
	}

	static void corners(const Mesh &mesh, const std::string &path, double shrinkSubdomain, double shringCluster)
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

	static void dirichlet(const Mesh &mesh, const std::string &path, double shrinkSubdomain, double shringCluster)
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

		outputx.store(dx, shrinkSubdomain, shringCluster);
		outputy.store(dy, shrinkSubdomain, shringCluster);
		outputz.store(dz, shrinkSubdomain, shringCluster);
	}

protected:
	void coordinatesDisplacement(const std::vector<std::vector<double> > &displacement, size_t dofs);
};

}
}



#endif /* OUTPUT_VTK_FULL_VTK_H_ */
