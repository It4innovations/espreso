
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
		switch (configuration.linear_elasticity_2D.library) {
		case SOLVER_LIBRARY::ESPRESO:
			instance = new LinearInstance<LinearElasticity2D>(configuration.linear_elasticity_2D.espreso, mesh);
			break;
		case SOLVER_LIBRARY::HYPRE:
			instance = new HypreInstance<LinearElasticity2D>(configuration.linear_elasticity_2D.hypre, mesh);
			break;
		};
		break;
	case PHYSICS::LINEAR_ELASTICITY_3D:
		switch (configuration.linear_elasticity_3D.library) {
		case SOLVER_LIBRARY::ESPRESO:
			instance = new LinearInstance<Elasticity3D>(configuration.linear_elasticity_3D.espreso, mesh);
			break;
		case SOLVER_LIBRARY::HYPRE:
			instance = new HypreInstance<Elasticity3D>(configuration.linear_elasticity_3D.hypre, mesh);
			break;
		};
		break;
	case PHYSICS::ADVECTION_DIFFUSION_2D:
		switch (configuration.advection_diffusion_2D.library) {
		case SOLVER_LIBRARY::ESPRESO:
			instance = new LinearInstance<AdvectionDiffusion2D>(configuration.advection_diffusion_2D.espreso, mesh);
			break;
		case SOLVER_LIBRARY::HYPRE:
			instance = new HypreInstance<AdvectionDiffusion2D>(configuration.advection_diffusion_2D.hypre, mesh);
			break;
		};
		break;
	case PHYSICS::ADVECTION_DIFFUSION_3D:
		switch (configuration.advection_diffusion_3D.library) {
		case SOLVER_LIBRARY::ESPRESO:
			instance = new LinearInstance<AdvectionDiffusion3D>(configuration.advection_diffusion_3D.espreso, mesh);
			break;
		case SOLVER_LIBRARY::HYPRE:
			instance = new HypreInstance<AdvectionDiffusion3D>(configuration.advection_diffusion_3D.hypre, mesh);
			break;
		};
		break;
	default:
		ESINFO(GLOBAL_ERROR) << "Not implemented physics";
	}

}



