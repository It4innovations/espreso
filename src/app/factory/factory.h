
#ifndef APP_FACTORY_FACTORY_H_
#define APP_FACTORY_FACTORY_H_

#include "esmesh.h"
#include "esassembler.h"

namespace espreso {

enum class GeneratorShape {
	CUBE,
	SPHERE,
	PLANE,
	CUBES
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
	LINEAR_ELASTICITY_2D,
	LINEAR_ELASTICITY_3D,
	TRANSIENT_ELASTICITY_2D,
	TRANSIENT_ELASTICITY_3D,
	ADVECTION_DIFFUSION_2D,
	ADVECTION_DIFFUSION_3D,
	STOKES
};

struct Factory {

	Factory(const ArgsConfiguration &configuration);
	~Factory()
	{
		delete instance;
	}

	void solve(const std::string &outputFile);

	double norm() const;

	Instance *instance;
	Mesh mesh;

private:
	std::vector<std::vector<double> > _solution;

	void readParameters(const ArgsConfiguration &configuration);

	ElementType eType;
	GeneratorShape shape;
	PhysicsAssembler physics;
};

}



#endif /* APP_FACTORY_FACTORY_H_ */
