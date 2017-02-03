
#include "assembler.h"
#include "old_physics/elasticity2d/assembler.h"
#include "old_physics/elasticity3d/assembler.h"
#include "old_physics/linear/advectiondiffusion2d/assembler.h"
#include "old_physics/linear/advectiondiffusion3d/assembler.h"

#include "instance/linear/instance.h"
#include "instance/dynamics/instance.h"
#include "instance/hypre/instance.h"
#include "instance/nonlinear/ssnm/instance.h"

#include "../config/globalconfiguration.h"

#include "../assembler/instance.h"
#include "../assembler/solver/linear.h"
#include "../assembler/physics/advectiondiffusion2d.h"

#include <numeric>

using namespace espreso;

class Mesh;

void Assembler::compose(const GlobalConfiguration &configuration, Instance* &instance, Mesh &mesh)
{
	instance = NULL;

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
		case SOLVER_LIBRARY::ESPRESO: {
			instance = new LinearInstance<AdvectionDiffusion2D, AdvectionDiffusion2DConfiguration>(configuration.advection_diffusion_2D, configuration.output, mesh);
			if (configuration.advection_diffusion_2D.newassembler) {
				std::vector<NewInstance*> instances;
				std::vector<NewPhysics*> physics;
				std::vector<LinearSolver*> linearSolvers;
				store::ResultStore* store;

				instances.push_back(new NewInstance(mesh.parts()));
				physics.push_back(new NewAdvectionDiffusion2D(&mesh, instances.front(), configuration.advection_diffusion_2D));
				linearSolvers.push_back(new LinearSolver(configuration.advection_diffusion_2D.espreso, instance->physics(), instance->constraints()));
				store = new store::VTK(configuration.output, mesh, "results");

				Linear solver(&mesh, physics, instances, linearSolvers, store);
				solver.init();
				std::vector<std::vector<double> > solution;
				solver.solve(solution);
				solver.finalize();
				exit(0);
			}
		} break;
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



