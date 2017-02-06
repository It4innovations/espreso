
#include "factory.h"

#include "../../assembler/step.h"
#include "../../assembler/instance.h"

#include "../../assembler/solver/linear.h"
#include "../../assembler/solver/newtonrhapson.h"

#include "../../assembler/physics/advectiondiffusion2d.h"
#include "../../solver/generic/LinearSolver.h"
#include "../../output/vtk/vtk.h"

#include "../../input/loader.h"
#include "../../assembler/assembler.h"
#include "../../assembler/instance/instance.h"
#include "../../assembler/old_physics/assembler.h"
#include "../../configuration/globalconfiguration.h"
#include "../../mesh/structures/mesh.h"
#include "../../mesh/settings/evaluator.h"
#include "../../mesh/elements/element.h"

namespace espreso {

Factory::Factory(const GlobalConfiguration &configuration)
: store(NULL), instance(NULL), mesh(new Mesh()), _newAssembler(false)
{
	input::Loader::load(configuration, *mesh, configuration.env.MPIrank, configuration.env.MPIsize);
	Assembler::compose(configuration, instance, *mesh);

	if (configuration.physics == PHYSICS::ADVECTION_DIFFUSION_2D && configuration.advection_diffusion_2D.newassembler) {
		_newAssembler = true;
		_instances.push_back(new Instance(mesh->parts()));
		_physics.push_back(new NewAdvectionDiffusion2D(mesh, _instances.front(), configuration.advection_diffusion_2D));
		_linearSolvers.push_back(new LinearSolver(configuration.advection_diffusion_2D.espreso, instance->physics(), instance->constraints()));
		store = new store::VTK(configuration.output, *mesh, "results");

		_solvers.push_back(new NewtonRhapson(mesh, _physics, _instances, _linearSolvers, store));
		for (size_t i = 0; i < configuration.advection_diffusion_2D.physics_solver.load_steps; i++) {
			loadSteps.push_back(_solvers.back());
		}
		meshPreprocessing();
	}
}

void Factory::meshPreprocessing()
{
	for (size_t i = 0; i < _physics.size(); i++) {

		switch (_linearSolvers.front()->configuration.method) {
		case ESPRESO_METHOD::TOTAL_FETI:
			_physics[i]->prepareTotalFETI();
			break;
		case ESPRESO_METHOD::HYBRID_FETI:
			switch (_linearSolvers.back()->configuration.B0_type) {
			case B0_TYPE::CORNERS:
				_physics[i]->prepareHybridTotalFETIWithCorners();
				break;
			case B0_TYPE::KERNELS:
				_physics[i]->prepareHybridTotalFETIWithKernels();
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

	std::for_each(_solvers.begin(), _solvers.end(), [] (Solver* solver) { delete solver; });
	std::for_each(_physics.begin(), _physics.end(), [] (Physics* physics) { delete physics; });
	std::for_each(_instances.begin(), _instances.end(), [] (Instance* instance) { delete instance; });
	std::for_each(_linearSolvers.begin(), _linearSolvers.end(), [] (LinearSolver* linearSolver) { delete linearSolver; });
	delete store;
}

void Factory::solve()
{
	if (!_newAssembler) {
		instance->init();
		instance->solve(_solution);
		instance->finalize();
	} else {
		Step step;
		for (size_t loadStep = 0; loadStep < loadSteps.size(); loadStep++) {
			step.load = loadStep;
			loadSteps[loadStep]->run(step);
			_solution = _instances.front()->primalSolution;
		}
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


