
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
	Hexahedron8(const permoncube::Settings &settings): _settings(settings) {};

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


#endif /* PM_HEXAHEDRON8_H_ */
