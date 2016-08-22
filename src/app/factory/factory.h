
#ifndef APP_FACTORY_FACTORY_H_
#define APP_FACTORY_FACTORY_H_

#include "esmesh.h"
#include "esinput.h"
#include "esoutput.h"
#include "esassemblers.h"

namespace espreso {

enum class GeneratorShape {
	CUBE,
	SPHERE,
	PLANE,
	CUBES
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

	Factory(const Configuration &configuration);
	~Factory()
	{
		return _mesh;
	}

	void solve();
	void store(const std::string &file);

	double norm() const;

	Instance *instance;
	Mesh mesh;

private:
	AssemblerBase *_assembler;
	std::vector<std::vector<double> > _solution;

	Mesh *_mesh;
	Mesh *_surface;
};

}



#endif /* APP_FACTORY_FACTORY_H_ */
