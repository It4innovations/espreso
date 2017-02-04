
#include "factory.h"

#include "../../assembler/step.h"
#include "../../assembler/instance.h"

#include "../../assembler/solver/linear.h"
#include "../../assembler/physics/advectiondiffusion2d.h"
#include "../../solver/generic/LinearSolver.h"
#include "../../output/vtk/vtk.h"

#include "../../assembler/instance/instance.h"
#include "../../input/loader.h"
#include "../../assembler/assembler.h"
#include "../../assembler/old_physics/assembler.h"
#include "../../config/globalconfiguration.h"
#include "../../mesh/structures/mesh.h"
#include "../../mesh/settings/evaluator.h"
#include "../../mesh/elements/element.h"

namespace espreso {

Factory::Factory(const GlobalConfiguration &configuration)
: solver(NULL), store(NULL), instance(NULL), mesh(new Mesh()), newAssembler(false)
{
	input::Loader::load(configuration, *mesh, configuration.env.MPIrank, configuration.env.MPIsize);
	Assembler::compose(configuration, instance, *mesh);

	if (configuration.physics == PHYSICS::ADVECTION_DIFFUSION_2D && configuration.advection_diffusion_2D.newassembler) {
		newAssembler = true;
		instances.push_back(new NewInstance(mesh->parts()));
		physics.push_back(new NewAdvectionDiffusion2D(mesh, instances.front(), configuration.advection_diffusion_2D));
		linearSolvers.push_back(new LinearSolver(configuration.advection_diffusion_2D.espreso, instance->physics(), instance->constraints()));
		store = new store::VTK(configuration.output, *mesh, "results");

		solver = new Linear(mesh, physics, instances, linearSolvers, store);
		meshPreprocessing();
	}
}

void Factory::meshPreprocessing()
{
	for (size_t i = 0; i < physics.size(); i++) {

		switch (linearSolvers[i]->configuration.method) {
		case ESPRESO_METHOD::TOTAL_FETI:
			physics[i]->prepareTotalFETI();
			break;
		case ESPRESO_METHOD::HYBRID_FETI:
			switch (linearSolvers[i]->configuration.B0_type) {
			case B0_TYPE::CORNERS:
				physics[i]->prepareHybridTotalFETIWithCorners();
				break;
			case B0_TYPE::KERNELS:
				physics[i]->prepareHybridTotalFETIWithKernels();
				break;
			default:
				ESINFO(GLOBAL_ERROR) << "Unknown type of B0";
			}
			break;
		default:
			ESINFO(GLOBAL_ERROR) << "Unknown FETI method";
		}

	}
	store->storeGeometry();
}

Factory::~Factory()
{
	if (instance != NULL) {
		delete instance;
	}
	delete mesh;
}

void Factory::solve()
{
	if (!newAssembler) {
		instance->init();
		instance->solve(_solution);
		instance->finalize();
	} else {
		solver->run();
		_solution = instances.front()->primalSolution;
	}
}

void Factory::check(const Results &configuration)
{
	double epsilon = 1e-2;
	auto norm = [&] () {
		double n = 0, sum = 0;
		for (size_t i = 0; i < _solution.size(); i++) {
			for (size_t j = 0; j < _solution[i].size(); j++) {
				n += _solution[i][j] * _solution[i][j];
			}
		}

		MPI_Allreduce(&n, &sum, 1, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
		return sqrt(sum);
	};

	if (configuration.norm != 0) {
		double nn = norm();
		ESTEST(EVALUATION)
			<< (fabs(nn - configuration.norm) > 1e-2 && !environment->MPIrank ? TEST_FAILED : TEST_PASSED)
			<< (fabs(nn - configuration.norm) > epsilon && !environment->MPIrank ? TEST_FAILED : TEST_PASSED)
			<< "Norm of the solution " << nn << " is not " << configuration.norm << ".";
	}

	auto evaluateProperty = [&] (const std::string &value, Property property, size_t DOF) {
		if (instance->physics().pointDOFs[DOF] == property && value.size()) {
			CoordinatesEvaluator evaluator(value, mesh->coordinates());
			for (size_t p = 0; p < mesh->parts(); p++) {
				for (size_t n = 0; n < mesh->coordinates().localSize(p); n++) {
					eslocal index = mesh->coordinates().localToCluster(p)[n];
					ESTEST(EVALUATION)
						<< (fabs(evaluator.evaluate(mesh->coordinates()[index]) - _solution[p][mesh->nodes()[index]->DOFIndex(p, DOF)]) > epsilon ? TEST_FAILED : TEST_PASSED)
						<< "Incorrect " << property << " of the solution.";
				}
			}
		}
	};

	for (size_t DOF = 0; DOF < instance->physics().pointDOFs.size(); DOF++) {
		evaluateProperty(configuration.displacement_x, Property::DISPLACEMENT_X, DOF);
		evaluateProperty(configuration.displacement_y, Property::DISPLACEMENT_Y, DOF);
		evaluateProperty(configuration.displacement_z, Property::DISPLACEMENT_Z, DOF);
		evaluateProperty(configuration.temperature   , Property::TEMPERATURE   , DOF);
		evaluateProperty(configuration.pressure      , Property::PRESSURE      , DOF);
	};
}

double Factory::norm() const
{
	double n = 0, sum = 0;
	for (size_t i = 0; i < _solution.size(); i++) {
		for (size_t j = 0; j < _solution[i].size(); j++) {
			n += _solution[i][j] * _solution[i][j];
		}
	}

	MPI_Allreduce(&n, &sum, 1, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
	return sqrt(sum);
}

}


