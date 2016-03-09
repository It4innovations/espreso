
#ifndef OUTPUT_VTK_FULL_VTK_H_
#define OUTPUT_VTK_FULL_VTK_H_

#include "../vtk.h"

namespace esoutput {

class VTK_Full: public VTK {

public:
	VTK_Full(const mesh::Mesh &mesh, const std::string &path): VTK(mesh, path) { };

	static void mesh(const mesh::Mesh &mesh, const std::string &path, double shrinkSubdomain, double shringCluster)
	{
		VTK_Full output(mesh, path);
		output.store(shrinkSubdomain, shringCluster);
	}

	static void fixPoints(const mesh::Mesh &mesh, const std::string &path, double shrinkSubdomain, double shringCluster)
	{
		VTK_Full output(mesh, path);
		output.store(mesh.getFixPoints(), shrinkSubdomain, shringCluster);
	}

	static void corners(const mesh::Mesh &mesh, const std::string &path, double shrinkSubdomain, double shringCluster)
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

protected:
	void coordinatesDisplacement(const std::vector<std::vector<double> > &displacement, size_t dofs);
};


}



#endif /* OUTPUT_VTK_FULL_VTK_H_ */
