
#ifndef PM_HEXAHEDRON20_H_
#define PM_HEXAHEDRON20_H_

#define Hexahedron20Subelements 1
#define Hexahedron20Subnodes 1

#include "esmesh.h"
#include "../../settings.h"
#include "../../generator.h"
#include "../../utils.h"

namespace permoncube {

class Hexahedron20 {

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

	static void fillGlobalBoundaries(
			const permoncube::Settings &settings,
			mesh::Boundaries &boundaries);

	static void clear() { _projection.clear(); };

	static eslocal subnodes[3];
	static eslocal subelements;

private:
	static std::vector<eslocal> _projection;
};

}




#endif /* PM_HEXAHEDRON20_H_ */
