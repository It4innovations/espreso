
#ifndef OUTPUT_VTK_FULL_VTK_H_
#define OUTPUT_VTK_FULL_VTK_H_

#include "../vtk.h"

namespace espreso {
namespace output {

class VTK_Full: public VTK {

public:
	VTK_Full(const Mesh &mesh, const std::string &path): VTK(mesh, path) { };

	static void mesh(const Mesh &mesh, const std::string &path, double shrinkSubdomain, double shringCluster);
	static void fixPoints(const Mesh &mesh, const std::string &path, double shrinkSubdomain, double shringCluster);
	static void corners(const Mesh &mesh, const std::string &path, double shrinkSubdomain, double shringCluster);
	static void dirichlet(const Mesh &mesh, const std::string &path, double shrinkSubdomain, double shringCluster);
	static void averaging(const Mesh &mesh, const std::string &path, double shrinkSubdomain, double shringCluster);

protected:
	void coordinatesDisplacement(const std::vector<std::vector<double> > &displacement);
};

}
}



#endif /* OUTPUT_VTK_FULL_VTK_H_ */
