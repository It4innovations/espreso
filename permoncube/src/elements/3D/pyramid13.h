
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
		return 	(wire(x) && (!odd(y) || !odd(z))) ||
				(wire(y) && (!odd(x) || !odd(z))) ||
				(wire(z) && (!odd(y) || !odd(x))) ||
				(odd(x) && odd(y) && odd(z)) ||
				(mid(x) && mid(y) && mid(z));
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

	inline static bool wire(const esglobal &x)
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
