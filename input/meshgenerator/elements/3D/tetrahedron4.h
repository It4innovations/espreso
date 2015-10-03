
#ifndef INPUT_MESHGENERATOR_ELEMENTS_3D_TETRAHEDRON4_H_
#define INPUT_MESHGENERATOR_ELEMENTS_3D_TETRAHEDRON4_H_

#define Tetrahedron4Subelements 6
#define Tetrahedron4Subnodes 0

#include <vector>
#include <cstring>

#include "esmesh.h"
#include "../../../meshgenerator/cube/settings.h"
#include "../../../meshgenerator/cube/utils.h"

namespace esinput {

class Tetrahedron4 {

public:
	Tetrahedron4(const esinput::CubeSettings &settings);

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


#endif /* INPUT_MESHGENERATOR_ELEMENTS_3D_TETRAHEDRON4_H_ */
