
#include "assembler.h"

#include "instance/linear/instance.h"
#include "instance/dynamics/instance.h"
#include "instance/hypre/instance.h"
#include "instance/nonlinear/ssnm/instance.h"

#include "physics/linear/advectiondiffusion2d/assembler.h"
#include "physics/linear/advectiondiffusion3d/assembler.h"

#include "physics/elasticity2d/assembler.h"
#include "physics/elasticity3d/assembler.h"
#include "../config/globalconfiguration.h"

#include <numeric>

using namespace espreso;

class Mesh;

void Assembler::compose(const GlobalConfiguration &configuration, Instance* &instance, Mesh &mesh)
{
	switch (configuration.physics) {
	case PHYSICS::LINEAR_ELASTICITY_2D:
		switch (configuration.linear_elasticity_2D.solver_library) {
		case SOLVER_LIBRARY::ESPRESO:
			instance = new LinearInstance<Elasticity2D, LinearElasticity2DConfiguration>(configuration.linear_elasticity_2D, configuration.output, mesh);
			break;
		case SOLVER_LIBRARY::HYPRE:
			instance = new HypreInstance<Elasticity2D, LinearElasticity2DConfiguration>(configuration.linear_elasticity_2D, configuration.output, mesh);
			break;
		};
		break;
	case PHYSICS::LINEAR_ELASTICITY_3D:
		switch (configuration.linear_elasticity_3D.solver_library) {
		case SOLVER_LIBRARY::ESPRESO:
			instance = new LinearInstance<Elasticity3D, LinearElasticity3DConfiguration>(configuration.linear_elasticity_3D, configuration.output, mesh);
			break;
		case SOLVER_LIBRARY::HYPRE:
			instance = new HypreInstance<Elasticity3D, LinearElasticity3DConfiguration>(configuration.linear_elasticity_3D, configuration.output, mesh);
			break;
		};
		break;
	case PHYSICS::ADVECTION_DIFFUSION_2D:
		switch (configuration.advection_diffusion_2D.solver_library) {
		case SOLVER_LIBRARY::ESPRESO:
			instance = new LinearInstance<AdvectionDiffusion2D, AdvectionDiffusion2DConfiguration>(configuration.advection_diffusion_2D, configuration.output, mesh);
			break;
		case SOLVER_LIBRARY::HYPRE:
			instance = new HypreInstance<AdvectionDiffusion2D, AdvectionDiffusion2DConfiguration>(configuration.advection_diffusion_2D, configuration.output, mesh);
			break;
		};
		break;
	case PHYSICS::ADVECTION_DIFFUSION_3D:
		switch (configuration.advection_diffusion_3D.solver_library) {
		case SOLVER_LIBRARY::ESPRESO:
			instance = new LinearInstance<AdvectionDiffusion3D, AdvectionDiffusion3DConfiguration>(configuration.advection_diffusion_3D, configuration.output, mesh);
			break;
		case SOLVER_LIBRARY::HYPRE:
			instance = new HypreInstance<AdvectionDiffusion3D, AdvectionDiffusion3DConfiguration>(configuration.advection_diffusion_3D, configuration.output, mesh);
			break;
		};
		break;
	default:
		ESINFO(GLOBAL_ERROR) << "Not implemented physics";
	}

}



