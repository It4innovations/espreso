
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
	Prisma6(const permoncube::Settings &settings): _settings(settings) {};

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



#endif /* PM_PRISMA6_H_ */
