
#ifndef OUTPUT_VTK_VTK_H_
#define OUTPUT_VTK_VTK_H_

#include "../store.h"

namespace esoutput {

class VTK: public ResultStore {

public:
	VTK(const mesh::Mesh &mesh, const std::string &path): ResultStore(mesh, path) { };

	void store(double shrinkSubdomain, double shringCluster);
	void store(std::vector<std::vector<double> > &displacement, size_t dofs, double shrinkSubdomain, double shringCluster);

protected:
	void head();
	virtual void coordinates(const mesh::Coordinates &coordinates, double shrinkSubdomain, double shringCluster);
	virtual void elements(const mesh::Mesh &mesh);
	virtual void coordinatesDisplacement(const std::vector<std::vector<double> > &displacement, size_t dofs) = 0;

	std::ofstream _vtk;
};


}


#endif /* OUTPUT_VTK_VTK_H_ */
