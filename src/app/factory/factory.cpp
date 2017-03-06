
#include "factory.h"

#include "../../assembler/step.h"
#include "../../assembler/instance.h"

#include "../../assembler/solver/linear.h"
#include "../../assembler/solver/newtonrhapson.h"

#include "../../assembler/solver/transientfirstorderimplicit.h"
#include "../../assembler/physics/advectiondiffusion2d.h"
#include "../../assembler/physics/shallowwater2d.h"
#include "../../solver/generic/LinearSolver.h"
#include "../../input/loader.h"
#include "../../assembler/assembler.h"
#include "../../assembler/instance/instance.h"
#include "../../assembler/old_physics/assembler.h"
#include "../../configuration/globalconfiguration.h"
#include "../../mesh/structures/mesh.h"
#include "../../mesh/settings/evaluator.h"
#include "../../mesh/elements/element.h"
#include "../../output/resultstorelist.h"

#include "../../output/resultstore/vtklegacy.h"
#include "../../output/resultstore/vtkxmlascii.h"
#include "../../output/resultstore/vtkxmlbinary.h"
#include "../../output/resultstore/catalyst.h"

namespace espreso {

Factory::Factory(const GlobalConfiguration &configuration)
: store(NULL), instance(NULL), mesh(new Mesh()), _newAssembler(false)
{
	input::Loader::load(configuration, *mesh, configuration.env.MPIrank, configuration.env.MPIsize);

	if (configuration.physics == PHYSICS::SHALLOW_WATER_2D) {
		_newAssembler = true;
		_instances.push_back(new Instance(mesh->parts(), mesh->neighbours()));
		_physics.push_back(new ShallowWater2D(mesh, _instances.front(), configuration.shallow_water_2D));
		_linearSolvers.push_back(new LinearSolver(_instances.front(), configuration.shallow_water_2D.espreso));

		loadSteps.push_back(new Linear(mesh, _physics.front(),  _linearSolvers.front(), store));
		meshPreprocessing(configuration.output);
		return;
	}

	Assembler::compose(configuration, instance, *mesh);

	store = new output::ResultStoreList(configuration.output);

	if (configuration.output.catalyst) {
		store->add(new output::Catalyst(configuration.output, mesh, "results"));
	}
	if (configuration.output.results || configuration.output.properties) {
		switch (configuration.output.format) {
		case OUTPUT_FORMAT::VTK_LEGACY:
			store->add(new output::VTKLegacy(configuration.output, mesh, "results"));
			break;
		case OUTPUT_FORMAT::VTK_XML_ASCII:
			store->add(new output::VTKXMLASCII(configuration.output, mesh, "results"));
			break;
		case OUTPUT_FORMAT::VTK_XML_BINARY:
			store->add(new output::VTKXMLBinary(configuration.output, mesh, "results"));
			break;
		default:
			ESINFO(GLOBAL_ERROR) << "ESPRESO internal error: add OUTPUT_FORMAT to factory.";
		}

	}

	if (configuration.physics == PHYSICS::ADVECTION_DIFFUSION_2D && configuration.advection_diffusion_2D.newassembler) {
		_newAssembler = true;
		_instances.push_back(new Instance(mesh->parts(), mesh->neighbours()));
		_physics.push_back(new NewAdvectionDiffusion2D(mesh, _instances.front(), configuration.advection_diffusion_2D));
		_linearSolvers.push_back(new LinearSolver(_instances.front(), configuration.advection_diffusion_2D.espreso));

		for (size_t i = 1; i <= configuration.advection_diffusion_2D.physics_solver.load_steps; i++) {
			auto it = configuration.advection_diffusion_2D.physics_solver.load_steps_settings.find(i);
			if (it == configuration.advection_diffusion_2D.physics_solver.load_steps_settings.end()) {
				loadSteps.push_back(new Linear(mesh, _physics.front(),  _linearSolvers.front(), store));
				break;
			}
			LoadStepSettings<AdvectionDiffusionNonLinearConvergence> *loadStepSettings = it->second;

			if (loadStepSettings->type == LoadStepSettingsBase::TYPE::STEADY_STATE) {
				if (loadStepSettings->mode == LoadStepSettingsBase::MODE::LINEAR) {
					loadSteps.push_back(new Linear(mesh, _physics.front(), _linearSolvers.front(), store));
				}
				if (loadStepSettings->mode == LoadStepSettingsBase::MODE::NONLINEAR) {
					switch (loadStepSettings->nonlinear_solver.method) {
					case NonLinearSolverBase::METHOD::NEWTON_RHAPSON:
						loadSteps.push_back(new NewtonRhapson(mesh, _physics.front(), _linearSolvers.front(), store, loadStepSettings->nonlinear_solver));
						break;
					case NonLinearSolverBase::METHOD::MODIFIED_NEWTON_RHAPSON:
						loadSteps.push_back(new NewtonRhapson(mesh, _physics.front(), _linearSolvers.front(), store, loadStepSettings->nonlinear_solver));
						break;
					default:
						ESINFO(GLOBAL_ERROR) << "Not implemented non-linear solver method";
					}
				}
			}
			if (loadStepSettings->type == LoadStepSettingsBase::TYPE::TRANSIENT) {
				if (loadStepSettings->mode == LoadStepSettingsBase::MODE::LINEAR) {
					switch (loadStepSettings->transient_solver.method) {
					case TransientSolver::METHOD::CRANK_NICOLSON:
					case TransientSolver::METHOD::GALERKIN:
					case TransientSolver::METHOD::BACKWARD_DIFF:
						loadSteps.push_back(new TransientFirstOrderImplicit(mesh, _physics.front(), _linearSolvers.front(), store, loadStepSettings->transient_solver, loadStepSettings->duration_time));
						break;
					case TransientSolver::METHOD::FORWARD_DIFF:
						ESINFO(GLOBAL_ERROR) << "Not implemented transient solver method.";
						break;
					default:
						ESINFO(GLOBAL_ERROR) << "Not implemented transient solver method";
					}
				}
			}

			if (i != loadSteps.size()) {
				ESINFO(GLOBAL_ERROR) << "Not implemented Physics solver";
			}
		}
		meshPreprocessing(configuration.output);
	}
}

void Factory::meshPreprocessing(const OutputConfiguration &configuration)
{
	for (size_t i = 0; i < _physics.size(); i++) {

		switch (_linearSolvers.front()->configuration.method) {
		case ESPRESO_METHOD::TOTAL_FETI:
			_physics[i]->prepareTotalFETI();
			break;
		case ESPRESO_METHOD::HYBRID_FETI:
			if (_linearSolvers.back()->configuration.method == ESPRESO_METHOD::HYBRID_FETI && !mesh->isContinuous()) {
				ESINFO(GLOBAL_ERROR) << "Do not use HYBRID FETI for non-continuous clusters.";
			}
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
	if (Test::report(EXPENSIVE)) {
		mesh->checkRegions(mesh->nodes());
	}
	if (configuration.properties) {
		store->storeSettings(mesh->steps());
	}
}

void Factory::finalize()
{
	// Detele store while finalizing because of Catalyst
	delete store;
}

Factory::~Factory()
{
	if (instance != NULL) {
		delete instance;
	}
	delete mesh;

	std::for_each(loadSteps.begin(), loadSteps.end(), [] (Solver* solver) { delete solver; });
	std::for_each(_physics.begin(), _physics.end(), [] (Physics* physics) { delete physics; });
	std::for_each(_instances.begin(), _instances.end(), [] (Instance* instance) { delete instance; });
	std::for_each(_linearSolvers.begin(), _linearSolvers.end(), [] (LinearSolver* linearSolver) { delete linearSolver; });
}

void Factory::solve()
{
	if (!_newAssembler) {
		instance->init();
		instance->solve(_solution);
		instance->finalize();
	} else {
		Step step;
		Logging::step = &step;
		for (size_t loadStep = 0; loadStep < loadSteps.size(); loadStep++) {
			step.step = loadStep;
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
		double nn;
		if (!_newAssembler) {
			nn = norm();
		} else {
			nn = sqrt(_physics.front()->sumSquares(_solution, Physics::SumOperation::AVERAGE));
		}
		ESTEST(EVALUATION)
			<< (fabs((nn - configuration.norm) / nn) > epsilon && !environment->MPIrank ? TEST_FAILED : TEST_PASSED)
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


