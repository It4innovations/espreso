
#ifndef PM_PRISMA15_H_
#define PM_PRISMA15_H_

#define Prisma15Subelements 2
#define Prisma15Subnodes 1

#include "esmesh.h"
#include "prisma6.h"
#include "../../settings.h"
#include "../../utils.h"

namespace permoncube {

class Prisma15 {

public:
	Prisma15(const permoncube::Settings &settings);

	void addElements(mesh::Mesh &mesh, const eslocal indices[]);
	static eslocal clusterNodesCount(const permoncube::Settings &settings);
	static esglobal globalNodesCount(const permoncube::Settings &settings);

	inline bool addPoint(const esglobal &x, const esglobal &y, const esglobal &z)
	{
		return !(odd(z) && (odd(y) || odd(x)));
	}

	inline eslocal projectPoint(const eslocal &index)
	{
		return _projection[index];
	}

	inline esglobal offset_x(esglobal x, esglobal y, esglobal z)
	{
		if (z % 2 == 0) {
			return x;
		} else {
			return x / 2;
		}
	}

	inline esglobal offset_y(esglobal y, esglobal z)
	{
		if (z % 2 == 0) {
			return _g3Nodes[0] * y;
		} else {
			return _g2Nodes[0] * ((y + 1) / 2);
		}
	}

	inline esglobal offset_z(esglobal z)
	{
		return _g3Nodes[0] * _g3Nodes[1] * ((z + 1) / 2) + _g2Nodes[0] * _g2Nodes[1] * (z / 2);
	}

	static eslocal subnodes[3];
	static eslocal subelements;

private:
	const permoncube::Settings &_settings;

	esglobal _g2Nodes[3];
	esglobal _g3Nodes[3];

	std::vector<eslocal> _projection;

	inline static bool odd(const esglobal &x)
	{
		return x % 2 == 1;
	}
};

}

#endif /* PM_PRISMA15_H_ */
