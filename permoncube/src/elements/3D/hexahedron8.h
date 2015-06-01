
#ifndef PM_HEXAHEDRON8_H_
#define PM_HEXAHEDRON8_H_

#define Hexahedron8Subelements 1
#define Hexahedron8Subnodes 0

#include "esmesh.h"
#include "../../settings.h"
#include "../../generator.h"
#include "../../utils.h"

namespace permoncube {

class Hexahedron8 {

public:
	static void addElements(mesh::Mesh &mesh, const esint indices[]);
	static void addCoordinates(mesh::Mesh &mesh, const permoncube::Settings &settings, const size_t cluster[]);
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


#endif /* PM_HEXAHEDRON8_H_ */
