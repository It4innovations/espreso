
#ifndef PM_TETRAHEDRON10_H_
#define PM_TETRAHEDRON10_H_

#define Tetrahedron10Subelements 6
#define Tetrahedron10Subnodes 1

#include "esmesh.h"
#include "../../settings.h"
#include "../../generator.h"
#include "../../utils.h"

namespace permoncube {

class Tetrahedron10 {

public:
	static void addElements(mesh::Mesh &mesh, const eslocal indices[]);
	static void addCoordinates(mesh::Mesh &mesh, const permoncube::Settings &settings, const size_t cluster[]);
	static void fixZeroPlanes(
			const permoncube::Settings &settings,
			std::map<eslocal, double> &dirichlet_x,
			std::map<eslocal, double> &dirichlet_y,
			std::map<eslocal, double> &dirichlet_z,
			const size_t cluster[]);
	static void fixBottom(
			const permoncube::Settings &settings,
			std::map<eslocal, double> &dirichlet_x,
			std::map<eslocal, double> &dirichlet_y,
			std::map<eslocal, double> &dirichlet_z,
			const size_t cluster[]);

	static void clear() { };

	static eslocal subnodes[3];
	static eslocal subelements;
};

}




#endif /* PM_TETRAHEDRON10_H_ */
