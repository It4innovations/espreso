
#ifndef PM_PYRAMID13_H_
#define PM_PYRAMID13_H_

#define Pyramid13Subelements 6
#define Pyramid13Subnodes 3

#include "esmesh.h"
#include "../../settings.h"
#include "../../generator.h"
#include "../../utils.h"

namespace permoncube {

class Pyramid13 {

public:
	Pyramid13(const permoncube::Settings &settings);

	void addElements(mesh::Mesh &mesh, const eslocal indices[]);
	static eslocal clusterNodesCount(const permoncube::Settings &settings);
	static esglobal globalNodesCount(const permoncube::Settings &settings);

	inline bool addPoint(const esglobal &x, const esglobal &y, const esglobal &z)
	{
		return	(face(z) && face(y) && face(x)) ||
				(face(z) && mid(y) && face(x)) ||
				(face(z) && face(y) && mid(x)) ||
				(mid(z) && face(y) && face(x)) ||
				(odd(x) && odd(y) && odd(z)) ||
				(mid(x) && mid(y) && mid(z));
	}

	inline eslocal projectPoint(const eslocal &index)
	{
		return _projection[index];
	}

	inline esglobal offset_x(esglobal x, esglobal y, esglobal z)
	{
		if (face(z)) {
			return face(y) ? x / 2 : x / 4;
		}
		if (mid(z)) {
			return x / 4;
		}
		return x / 2;
	}

	inline esglobal offset_y(esglobal y, esglobal z)
	{
		if (face(z)) {
			return _g3Nodes[0] * ((y + 3) / 4) + _g2Nodes[0] * ((y + 1) / 4);
		}
		if (mid(z)) {
			return _g2Nodes[0] * ((y + 3) / 4) + (_g2Nodes[0] - 1) * ((y + 1) / 4);
		}
		return _g2Nodes[0] * (y / 2);
	}

	inline esglobal offset_z(esglobal z)
	{
		return
			faceNodes * ((z + 3) / 4) +
			_g2Nodes[0] * _g2Nodes[1] * ((z + 1) / 2) +
			(_g2Nodes[0] * _g2Nodes[1] + (_g2Nodes[0] - 1) * (_g2Nodes[1] - 1)) * ((z + 1) / 3);
	}

	static eslocal subnodes[3];
	static eslocal subelements;

private:
	const permoncube::Settings &_settings;

	std::vector<eslocal> _projection;

	esglobal _g3Nodes[3];
	esglobal _g2Nodes[3];
	esglobal faceNodes;

	inline static bool odd(const esglobal &x)
	{
		return x % 2 == 1;
	}

	inline static bool face(const esglobal &x)
	{
		return x % 4 == 0;
	}

	inline static bool mid(const esglobal &x)
	{
		return (x + 2) % 4 == 0;
	}
};

}
#endif /* PM_PYRAMID13_H_ */
