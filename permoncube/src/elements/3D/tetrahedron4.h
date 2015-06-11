
#ifndef PM_TETRAHEDRON4_H_
#define PM_TETRAHEDRON4_H_

#define Tetrahedron4Subelements 6
#define Tetrahedron4Subnodes 0

#include <vector>
#include <cstring>

#include "esmesh.h"
#include "../../settings.h"
#include "../../generator.h"
#include "../../utils.h"

namespace permoncube {

class Tetrahedron4 {

public:
	Tetrahedron4(const permoncube::Settings &settings): _settings(settings) {};

	void addElements(mesh::Mesh &mesh, const eslocal indices[]);
	static eslocal clusterNodesCount(const permoncube::Settings &settings);
	static esglobal globalNodesCount(const permoncube::Settings &settings);

	inline bool addPoint(const esglobal &x, const esglobal &y, const esglobal &z)
	{
		return true;
	}

	inline eslocal projectPoint(const eslocal &index)
	{
		return index;
	}

	static eslocal subnodes[3];
	static eslocal subelements;

private:
	const permoncube::Settings &_settings;
};

}


#endif /* PM_TETRAHEDRON4_H_ */
