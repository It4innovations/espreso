
#ifndef INPUT_MESHGENERATOR_ELEMENTS_3D_PYRAMID5_H_
#define INPUT_MESHGENERATOR_ELEMENTS_3D_PYRAMID5_H_

#define Pyramid5Subelements 6
#define Pyramid5Subnodes 1

#include "esmesh.h"
#include "../../../meshgenerator/elements/3D/hexahedron8.h"
#include "../../cube/settings.h"
#include "../../cube/utils.h"

namespace esinput {

class Pyramid5 {

public:
	Pyramid5(const esinput::CubeSettings &settings);

	void addElements(mesh::Mesh &mesh, const eslocal indices[]);
	static eslocal clusterNodesCount(const esinput::CubeSettings &settings);
	static esglobal globalNodesCount(const esinput::CubeSettings &settings);

	inline bool addPoint(const esglobal &x, const esglobal &y, const esglobal &z)
	{
		return (odd(x) && odd(y) && odd(z)) || (!odd(x) && !odd(y) && !odd(z));
	}

	inline eslocal projectPoint(const eslocal &index)
	{
		return _projection[index];
	}

	inline esglobal offset_x(esglobal x, esglobal y, esglobal z)
	{
		return x / 2;
	}

	inline esglobal offset_y(esglobal y, esglobal z)
	{
		if (z % 2 == 0) {
			return _gNodes[0] * ((y + 1) / 2);
		} else {
			return (_gNodes[0] - 1) * (y / 2);
		}
	}

	inline esglobal offset_z(esglobal z)
	{
		return _gNodes[0] * _gNodes[1] * ((z + 1) / 2) + (_gNodes[0] - 1) * (_gNodes[1] - 1) * (z / 2);
	}

	static eslocal subnodes[3];
	static eslocal subelements;

private:
	const esinput::CubeSettings &_settings;

	esglobal _gNodes[3];

	std::vector<eslocal> _projection;

	inline static bool odd(const esglobal &x)
	{
		return x % 2 == 1;
	}
};

}


#endif /* INPUT_MESHGENERATOR_ELEMENTS_3D_PYRAMID5_H_ */
