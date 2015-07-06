
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
	Prisma6(const permoncube::Settings &settings);

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



#endif /* PM_PRISMA6_H_ */
