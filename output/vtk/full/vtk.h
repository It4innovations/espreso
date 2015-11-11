
#ifndef OUTPUT_VTK_FULL_VTK_H_
#define OUTPUT_VTK_FULL_VTK_H_

#include "../vtk.h"

namespace esoutput {

class VTK_Full: public VTK {

public:
	VTK_Full(const mesh::Mesh &mesh, const std::string &path): VTK(mesh, path) { };

protected:
	void coordinatesDisplacement(const std::vector<std::vector<double> > &displacement, size_t dofs);
};


}



#endif /* OUTPUT_VTK_FULL_VTK_H_ */
