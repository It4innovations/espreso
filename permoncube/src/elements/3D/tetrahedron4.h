
#ifndef PM_TETRAHEDRON4_H_
#define PM_TETRAHEDRON4_H_

#define Tetrahedron4Subelements 6
#define Tetrahedron4Subnodes 0

#include <vector>
#include <cstring>

#include "esmesh.h"
#include "element3D.h"
#include "../../settings.h"
#include "../../generator.h"
#include "../../utils.h"

namespace permoncube {

class Tetrahedron4 {

public:
	static void addElements(mesh::Mesh &mesh, const esint indices[]);
	static void addCoordinates(mesh::Mesh &mesh, const Settings &settings, const size_t cluster[]);
	static void fixZeroPlanes(
			const permoncube::Settings &settings,
			std::map<esint, double> &dirichlet_x,
			std::map<esint, double> &dirichlet_y,
			std::map<esint, double> &dirichlet_z,
			const size_t cluster[]);
	static void fixBottom(
			const permoncube::Settings &settings,
			std::map<esint, double> &dirichlet_x,
			std::map<esint, double> &dirichlet_y,
			std::map<esint, double> &dirichlet_z,
			const size_t cluster[]);

	static void clear() { };

	static esint subnodes[3];
	static esint subelements;
};

}


#endif /* PM_TETRAHEDRON4_H_ */
