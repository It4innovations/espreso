
#include "assembler.h"

#include "../config/description.h"
#include "instance/linear/instance.h"
#include "instance/dynamics/instance.h"
#include "instance/hypre/instance.h"
#include "instance/nonlinear/ssnm/instance.h"

#include "physics/linear/advectiondiffusion2d/assembler.h"
#include "physics/linear/advectiondiffusion3d/assembler.h"
#include "physics/linear/elasticity2d/assembler.h"

#include "physics/elasticity3d/assembler.h"

using namespace espreso;

class Mesh;

void Assembler::compose(const GlobalConfiguration &configuration, Instance* &instance, Mesh &mesh)
{
	switch (configuration.physics) {
	case PHYSICS::LINEAR_ELASTICITY_2D:
		instance = new LinearInstance<LinearElasticity2D>(mesh);
		break;
	case PHYSICS::LINEAR_ELASTICITY_3D:
		instance = new LinearInstance<Elasticity3D>(mesh);
		break;
	case PHYSICS::ADVECTION_DIFFUSION_2D:
		instance = new LinearInstance<AdvectionDiffusion2D>(mesh);
		break;
	case PHYSICS::ADVECTION_DIFFUSION_3D:
		instance = new LinearInstance<AdvectionDiffusion3D>(mesh);
		break;
	default:
		ESINFO(GLOBAL_ERROR) << "Not implemented physics";
	}

}



