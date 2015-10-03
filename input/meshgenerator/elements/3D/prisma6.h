
#ifndef INPUT_MESHGENERATOR_ELEMENTS_3D_PRISMA6_H_
#define INPUT_MESHGENERATOR_ELEMENTS_3D_PRISMA6_H_

#define Prisma6Subelements 2
#define Prisma6Subnodes 0

#include "esmesh.h"
#include "../../../meshgenerator/cube/settings.h"
#include "../../../meshgenerator/cube/utils.h"

namespace esinput {

class Prisma6 {

public:
	Prisma6(const esinput::CubeSettings &settings);

	void addElements(mesh::Mesh &mesh, const eslocal indices[]);
	static eslocal clusterNodesCount(const esinput::CubeSettings &settings);
	static esglobal globalNodesCount(const esinput::CubeSettings &settings);

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
	const esinput::CubeSettings &_settings;

	esglobal _gNodes[3];
};

}



#endif /* INPUT_MESHGENERATOR_ELEMENTS_3D_PRISMA6_H_ */
