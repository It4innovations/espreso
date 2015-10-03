
#ifndef INPUT_MESHGENERATOR_ELEMENTS_3D_HEXAHEDRON20_H_
#define INPUT_MESHGENERATOR_ELEMENTS_3D_HEXAHEDRON20_H_

#define Hexahedron20Subelements 1
#define Hexahedron20Subnodes 1

#include "esmesh.h"
#include "../../../meshgenerator/elements/3D/hexahedron8.h"
#include "../../cube/settings.h"
#include "../../cube/utils.h"

namespace esinput {

class Hexahedron20 {

public:
	Hexahedron20(const esinput::CubeSettings &settings);

	void addElements(mesh::Mesh &mesh, const eslocal indices[]);
	static eslocal clusterNodesCount(const esinput::CubeSettings &settings);
	static esglobal globalNodesCount(const esinput::CubeSettings &settings);

	inline bool addPoint(const esglobal &x, const esglobal &y, const esglobal &z)
	{
		return !((odd(x) && odd(y)) || (odd(y) && odd(z)) || (odd(x) && odd(z)));
	}

	inline eslocal projectPoint(const eslocal &index)
	{
		return _projection[index];
	}

	inline esglobal offset_x(esglobal x, esglobal y, esglobal z)
	{
		if (z % 2 == 0 && y % 2 == 0) {
			return x;
		} else {
			return x / 2;
		}
	}

	inline esglobal offset_y(esglobal y, esglobal z)
	{
		if (z % 2 == 0) {
			return _g3Nodes[0] * ((y + 1) / 2) + _g2Nodes[0] * (y / 2);
		} else {
			return _g2Nodes[0] * ((y + 1) / 2);
		}
	}

	inline esglobal offset_z(esglobal z)
	{
		return  faceNodes * ((z + 1) / 2) + _g2Nodes[0] * _g2Nodes[1] * (z / 2);
	}

	static eslocal subnodes[3];
	static eslocal subelements;

private:
	const esinput::CubeSettings &_settings;

	esglobal _g2Nodes[3];
	esglobal _g3Nodes[3];
	esglobal faceNodes;

	std::vector<eslocal> _projection;

	inline static bool odd(const esglobal &x)
	{
		return x % 2 == 1;
	}
};

}


#endif /* INPUT_MESHGENERATOR_ELEMENTS_3D_HEXAHEDRON20_H_ */
