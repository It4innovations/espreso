
#ifndef PM_PRISMA6_H_
#define PM_PRISMA6_H_

#define Prisma6Subelements 2
#define Prisma6Subnodes 0

#include "esmesh.h"
#include "../../settings.h"
#include "../../generator.h"
#include "../../utils.h"

namespace permoncube {

class Prisma6 {

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

	static void clear() { };

	static eslocal subnodes[3];
	static eslocal subelements;
};

}



#endif /* PM_PRISMA6_H_ */
