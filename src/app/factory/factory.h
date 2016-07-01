
#ifndef APP_FACTORY_FACTORY_H_
#define APP_FACTORY_FACTORY_H_

#include "esmesh.h"
#include "esassembler.h"

namespace espreso {

enum class GeneratorShape {
	CUBE,
	SPHERE,
	PLANE
};

enum class ElementType {
	HEXA8,
	HEXA20,
	TETRA4,
	TETRA10,
	PRISMA6,
	PRISMA15,
	PYRAMID5,
	PYRAMID13,

	SQUARE4,
	SQUARE8,
	TRIANGLE3,
	TRIANGLE6
};

enum class PhysicsAssembler {
	LINEAR_ELASTICITY,
	TEMPERATURE,
	TRANSIENT_ELASTICITY,
	ADVECTION_DIFFUSION
};

struct Factory {

	Factory(const Configuration &configuration);
	~Factory()
	{
		delete instance;
	}

	void solve(const std::string &outputFile);

	Instance *instance;
	Mesh mesh;

private:
	std::vector<std::vector<double> > _solution;

	void readParameters(const Configuration &configuration);

	ElementType eType;
	GeneratorShape shape;
	PhysicsAssembler physics;
};

}



#endif /* APP_FACTORY_FACTORY_H_ */
