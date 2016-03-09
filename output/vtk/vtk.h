
#ifndef OUTPUT_VTK_VTK_H_
#define OUTPUT_VTK_VTK_H_

#include "../store.h"

namespace esoutput {

class VTK: public ResultStore {

public:
	VTK(const mesh::Mesh &mesh, const std::string &path);

	virtual void store(const std::vector<std::vector<eslocal> > &points, double shrinkSubdomain, double shringCluster);

	void store(double shrinkSubdomain, double shringCluster);
	void store(std::vector<std::vector<double> > &displacement, size_t dofs, double shrinkSubdomain, double shringCluster);

protected:
	static mesh::Point shrink(const mesh::Point &p,
			const mesh::Point &subdomainCenter, double subdomainShrinkRatio,
			const mesh::Point &clusterCenter, double clusterShrinkRatio);

	void head();
	virtual void coordinates(const mesh::Coordinates &coordinates, double shrinkSubdomain, double shringCluster);
	virtual void coordinates(const mesh::Coordinates &coordinates, const std::vector<std::vector<eslocal> > &points, double shrinkSubdomain, double shringCluster);
	virtual void elements(const mesh::Mesh &mesh);
	virtual void points(const std::vector<std::vector<eslocal> > &points);
	virtual void coordinatesDisplacement(const std::vector<std::vector<double> > &displacement, size_t dofs) = 0;

	std::ofstream _vtk;

	std::vector<mesh::Point> _subdomainsCenter;
	mesh::Point _clusterCenter;
};


}


#endif /* OUTPUT_VTK_VTK_H_ */
