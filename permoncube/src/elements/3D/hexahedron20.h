
#ifndef PM_HEXAHEDRON20_H_
#define PM_HEXAHEDRON20_H_

#define Hexahedron20Subelements 1
#define Hexahedron20Subnodes 1

#include "esmesh.h"
#include "../../settings.h"
#include "../../generator.h"
#include "../../utils.h"

namespace permoncube {

class Hexahedron20 {

public:
	Hexahedron20(const permoncube::Settings &settings);

	void addElements(mesh::Mesh &mesh, const eslocal indices[]);
	static eslocal clusterNodesCount(const permoncube::Settings &settings);
	static esglobal globalNodesCount(const permoncube::Settings &settings);

	inline bool addPoint(const esglobal &x, const esglobal &y, const esglobal &z)
	{
		return !((odd(x) && odd(y)) || (odd(y) && odd(z)) || (odd(x) && odd(z)));
	}

	inline eslocal projectPoint(const eslocal &index)
	{
		return _projection[index];
	}

	static eslocal subnodes[3];
	static eslocal subelements;

private:
	const permoncube::Settings &_settings;

	std::vector<eslocal> _projection;

	inline static bool odd(const esglobal &x)
	{
		return x % 2 == 1;
	}
};

}




#endif /* PM_HEXAHEDRON20_H_ */
