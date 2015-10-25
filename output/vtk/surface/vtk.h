
#ifndef OUTPUT_VTK_SURFACE_VTK_H_
#define OUTPUT_VTK_SURFACE_VTK_H_

#include "../vtk.h"

namespace esoutput {

class VTK_Surface: public VTK {

public:
	VTK_Surface(const mesh::Mesh &mesh, const std::string &path): _full(mesh), _surface(mesh), VTK(_surface, path)  { };

protected:
	void coordinatesDisplacement(const std::vector<std::vector<double> > &displacement);

	const mesh::Mesh &_full;
	mesh::SurfaceMesh _surface;
};


}




#endif /* OUTPUT_VTK_SURFACE_VTK_H_ */
