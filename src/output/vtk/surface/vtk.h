
#ifndef OUTPUT_VTK_SURFACE_VTK_H_
#define OUTPUT_VTK_SURFACE_VTK_H_

#include "../vtk.h"

namespace espreso {
namespace output {

class VTK_Surface: public VTK {

public:
	VTK_Surface(const Mesh &mesh, const std::string &path): _full(mesh), VTK(_surface, path)
	{
		mesh.getSurface(_surface);
	};

protected:
	void coordinatesDisplacement(const std::vector<std::vector<double> > &displacement, size_t dofs);

	const Mesh &_full;
	Mesh _surface;
};

}
}




#endif /* OUTPUT_VTK_SURFACE_VTK_H_ */
