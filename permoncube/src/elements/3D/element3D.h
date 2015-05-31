
#ifndef PM_ELEMENT3D_H_
#define PM_ELEMENT3D_H_

#include "esmesh.h"
#include "../../settings.h"
#include "../../utils.h"

namespace permoncube {

template<class TElement>
class Element3D {

public:
	static void addFullCoordinates(mesh::Mesh &mesh, const Settings &settings, const size_t cluster[], std::vector<idx_t> &mapping);

	static void fixFullBottom(
			const permoncube::Settings &settings,
			std::map<int, double> &dirichlet_x,
			std::map<int, double> &dirichlet_y,
			std::map<int, double> &dirichlet_z,
			const size_t cluster[],
			std::vector<idx_t> &mapping);

	static void fixFullZeroPlanes(
			const permoncube::Settings &settings,
			std::map<int, double> &dirichlet_x,
			std::map<int, double> &dirichlet_y,
			std::map<int, double> &dirichlet_z,
			const size_t cluster[],
			std::vector<idx_t> &mapping);

};
}

#include "element3D.hpp"


#endif /* PM_ELEMENT3D_H_ */
