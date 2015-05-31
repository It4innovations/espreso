
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
	static void addElements(mesh::Mesh &mesh, const idx_t indices[]);
	static void addCoordinates(mesh::Mesh &mesh, const permoncube::Settings &settings, const size_t cluster[]);
	static void fixZeroPlanes(
			const permoncube::Settings &settings,
			std::map<int, double> &dirichlet_x,
			std::map<int, double> &dirichlet_y,
			std::map<int, double> &dirichlet_z,
			const size_t cluster[]);
	static void fixBottom(
			const permoncube::Settings &settings,
			std::map<int, double> &dirichlet_x,
			std::map<int, double> &dirichlet_y,
			std::map<int, double> &dirichlet_z,
			const size_t cluster[]);

	static void clear();

	static size_t subnodes[3];
	static size_t subelements;

private:
	static std::vector<idx_t> _coordinateMapping;
};

}


#endif /* PM_HEXAHEDRON8_H_ */
