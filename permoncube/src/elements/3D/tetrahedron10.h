
#ifndef PM_TETRAHEDRON10_H_
#define PM_TETRAHEDRON10_H_

#define Tetrahedron10Subelements 6
#define Tetrahedron10Subnodes 1

#include "esmesh.h"
#include "../../settings.h"
#include "../../generator.h"

namespace permoncube {

class Tetrahedron10 {

public:
	static void addElements(mesh::Mesh &mesh, const idx_t indices[]);
	static void addCoordinates(mesh::Mesh &mesh, const permoncube::Settings &settings, const size_t cluster[]);
	static void clear();

	static size_t subnodes[3];
	static size_t subelements;

private:
	static std::vector<idx_t> _coordinateMapping;
};

}




#endif /* PM_TETRAHEDRON10_H_ */
