
#ifndef OUTPUT_VTK_VTK_H_
#define OUTPUT_VTK_VTK_H_

#include "../store.h"

namespace espreso {
namespace output {

class VTK: public ResultStore {

public:
	VTK(const Mesh &mesh, const std::string &path);

	virtual void store(const std::vector<std::vector<eslocal> > &points, double shrinkSubdomain, double shringCluster);

	void store(double shrinkSubdomain, double shringCluster);
	void store(std::vector<std::vector<double> > &displacement, size_t dofs, double shrinkSubdomain, double shringCluster);

protected:
	static Point shrink(const Point &p,
			const Point &subdomainCenter, double subdomainShrinkRatio,
			const Point &clusterCenter, double clusterShrinkRatio);

	void head();
	virtual void coordinates(const Coordinates &coordinates, double shrinkSubdomain, double shringCluster);
	virtual void coordinates(const Coordinates &coordinates, const std::vector<std::vector<eslocal> > &points, double shrinkSubdomain, double shringCluster);
	virtual void elements(const Mesh &mesh);
	virtual void points(const std::vector<std::vector<eslocal> > &points);
	virtual void coordinatesDisplacement(const std::vector<std::vector<double> > &displacement, size_t dofs) = 0;

	std::ofstream _vtk;

	std::vector<Point> _subdomainsCenter;
	Point _clusterCenter;
};


}
}


#endif /* OUTPUT_VTK_VTK_H_ */
