
#ifndef PM_ELEMENT3D_H_
#define PM_ELEMENT3D_H_

#include "esmesh.h"
#include "../../settings.h"
#include "../../utils.h"

namespace permoncube {

template<class TElement>
class Element3D {

public:
	static void addFullCoordinates(mesh::Mesh &mesh, const Settings &settings, const size_t cluster[]);

	static void fixFullBottom(
			const permoncube::Settings &settings,
			std::map<eslocal, double> &dirichlet_x,
			std::map<eslocal, double> &dirichlet_y,
			std::map<eslocal, double> &dirichlet_z,
			const size_t cluster[]);

	static void fixFullZeroPlanes(
			const permoncube::Settings &settings,
			std::map<eslocal, double> &dirichlet_x,
			std::map<eslocal, double> &dirichlet_y,
			std::map<eslocal, double> &dirichlet_z,
			const size_t cluster[]);

};
}

#include "element3D.hpp"


#endif /* PM_ELEMENT3D_H_ */
