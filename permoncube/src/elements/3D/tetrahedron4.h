
#ifndef PM_TETRAHEDRON4_H_
#define PM_TETRAHEDRON4_H_

#define Tetrahedron4Subelements 6
#define Tetrahedron4Subnodes 0

#include <vector>
#include <cstring>

#include "esmesh.h"
#include "../../settings.h"
#include "../../generator.h"

namespace permoncube {

class Tetrahedron4 {

public:
	static void addElements(mesh::Mesh &mesh, const idx_t indices[]);
	static void addCoordinates(mesh::Mesh &mesh, const Settings &settings, const size_t cluster[]);
	static void fixZeroPlanes(
			const permoncube::Settings &settings,
			std::map<int, double> &dirichlet_x,
			std::map<int, double> &dirichlet_y,
			std::map<int, double> &dirichlet_z);
	static void fixBottom(
			const permoncube::Settings &settings,
			std::map<int, double> &dirichlet_x,
			std::map<int, double> &dirichlet_y,
			std::map<int, double> &dirichlet_z);

	static void clear();

	static size_t subnodes[3];
	static size_t subelements;

private:
	static std::vector<idx_t> _coordinateMapping;
};

}


#endif /* PM_TETRAHEDRON4_H_ */
