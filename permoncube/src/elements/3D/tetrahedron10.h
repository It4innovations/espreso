
#ifndef PM_TETRAHEDRON10_H_
#define PM_TETRAHEDRON10_H_

#define Tetrahedron10Subelements 6
#define Tetrahedron10Subnodes 1

#include "esmesh.h"
#include "../../settings.h"
#include "../../generator.h"
#include "../../utils.h"

namespace permoncube {

class Tetrahedron10 {

public:
	Tetrahedron10(const permoncube::Settings &settings): _settings(settings) {};

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

	static eslocal subnodes[3];
	static eslocal subelements;

private:
	const permoncube::Settings &_settings;
};

}




#endif /* PM_TETRAHEDRON10_H_ */
