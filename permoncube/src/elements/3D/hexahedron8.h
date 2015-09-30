
#ifndef PM_HEXAHEDRON8_H_
#define PM_HEXAHEDRON8_H_

#define Hexahedron8Subelements 1
#define Hexahedron8Subnodes 0

#include "esmesh.h"
#include "../../settings.h"
#include "../../utils.h"

namespace permoncube {

class Hexahedron8 {

public:
	Hexahedron8(const permoncube::Settings &settings);

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

	inline esglobal offset_x(esglobal x, esglobal y, esglobal z)
	{
		return x;
	}

	inline esglobal offset_y(esglobal y, esglobal z)
	{
		return _gNodes[0] * y;
	}

	inline esglobal offset_z(esglobal z)
	{
		return _gNodes[0] * _gNodes[1] * z;
	}

	static eslocal subnodes[3];
	static eslocal subelements;

private:
	const permoncube::Settings &_settings;

	esglobal _gNodes[3];
};

}


#endif /* PM_HEXAHEDRON8_H_ */
