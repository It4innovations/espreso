
#include "factory.h"

#include "advectiondiffusionfactory.h"
#include "structuralmechanicsfactory.h"

#include "../../assembler/physicssolver/timestep/timestepsolver.h"
#include "../../assembler/physicssolver/loadstep/loadstepsolver.h"
#include "../../assembler/physicssolver/assembler.h"
#include "../../assembler/physics/physics.h"
#include "../../assembler/step.h"
#include "../../assembler/instance.h"
#include "../../input/loader.h"

#include "../../configuration/globalconfiguration.h"
#include "../../mesh/structures/mesh.h"

#include "../../output/resultstorelist.h"
#include "../../output/resultstore/asyncstore.h"
#include "../../output/resultstore/vtklegacy.h"
#include "../../output/resultstore/vtkxmlascii.h"
#include "../../output/resultstore/vtkxmlbinary.h"
#include "../../output/resultstore/catalyst.h"
#include "../../output/monitoring/monitoring.h"

#include "../../solver/generic/FETISolver.h"


namespace espreso {

Factory::Factory(const GlobalConfiguration &configuration)
: _mesh(new Mesh()), _storeList(new output::ResultStoreList(configuration.output))
{
	initAsync(configuration.output);

	_dispatcher.init();
	environment->MPICommunicator = _dispatcher.commWorld();
	MPI_Comm_rank(environment->MPICommunicator, &environment->MPIrank);
	MPI_Comm_size(environment->MPICommunicator, &environment->MPIsize);

	if (!_dispatcher.dispatch()) {
		return;
	}

	input::Loader::load(configuration, *_mesh, configuration.env.MPIrank, configuration.env.MPIsize);

	loadPhysics(configuration);
	setOutput(configuration.output);
}

void Factory::solve()
{
	if (_dispatcher.isExecutor()) {
		return;
	}

	Step step;
	Logging::step = &step;

	for (step.step = 0; step.step < _loadSteps.size(); step.step++) {
		_loadSteps[step.step]->run(step);
	}
}

void Factory::initAsync(const OutputConfiguration &configuration)
{
	_asyncStore = NULL;
	async::Config::setMode(async::SYNC);
	if (configuration.solution) {
		if (configuration.mode != OUTPUT_MODE::SYNC && (configuration.settings || configuration.FETI_data)) {
			ESINFO(ALWAYS) << Info::TextColor::YELLOW << "Storing of SETTINGS or FETI_DATA is implemented only for OUTPUT::MODE==SYNC. Hence, output is synchronized!";
		} else if (configuration.collected) {
			ESINFO(ALWAYS) << Info::TextColor::YELLOW << "Storing COLLECTED output is implemented only for OUTPUT::MODE==SYNC. Hence, output is synchronized!";
		} else {
			// Configure the asynchronous library
			switch (configuration.mode) {
			case OUTPUT_MODE::SYNC:
				async::Config::setMode(async::SYNC);
				break;
			case OUTPUT_MODE::THREAD:
				async::Config::setMode(async::THREAD);
				break;
			case OUTPUT_MODE::MPI:
				if (environment->MPIsize == 1) {
					ESINFO(GLOBAL_ERROR) << "Invalid number of MPI processes. OUTPUT::MODE==MPI required at least two MPI processes.";
				}
				if (configuration.output_node_group_size == 0) {
					ESINFO(GLOBAL_ERROR) << "OUTPUT::OUTPUT_NODE_GROUP_SIZE cannot be 0.";
				}
				async::Config::setMode(async::MPI);
				async::Config::setGroupSize(configuration.output_node_group_size);
				_dispatcher.setGroupSize(configuration.output_node_group_size);
				// async::Config::setUseAsyncCopy(true);
				break;
			}

			_asyncStore = new output::AsyncStore(configuration, _mesh);
		}
	}
}

void Factory::setOutput(const OutputConfiguration &configuration)
{
	if (configuration.catalyst) {
		_storeList->add(new output::Catalyst(configuration, _mesh));
	}
	if (configuration.solution || configuration.settings) {
		if (_asyncStore != NULL) {
			_asyncStore->init(_mesh);
			_storeList->add(_asyncStore);
		} else {
			switch (configuration.format) {
			case OUTPUT_FORMAT::VTK_LEGACY:
				_storeList->add(new output::VTKLegacy(configuration, _mesh));
				break;
			case OUTPUT_FORMAT::VTK_XML_ASCII:
				_storeList->add(new output::VTKXMLASCII(configuration, _mesh));
				break;
			case OUTPUT_FORMAT::VTK_XML_BINARY:
				_storeList->add(new output::VTKXMLBinary(configuration, _mesh));
				break;
			default:
				ESINFO(GLOBAL_ERROR) << "ESPRESO internal error: add OUTPUT_FORMAT to factory.";
			}
		}
	}
	if (configuration.monitoring.size()) {
		_storeList->add(new output::Monitoring(configuration, _mesh));
	}

	if (configuration.settings) {
		_storeList->storeSettings(_mesh->steps());
	}
}

FactoryLoader* Factory::createFactoryLoader(const GlobalConfiguration &configuration)
{
	switch (configuration.physics) {
	case PHYSICS::ADVECTION_DIFFUSION_2D:
		return new AdvectionDiffusionFactory(configuration.advection_diffusion_2D, _mesh);
	case PHYSICS::ADVECTION_DIFFUSION_3D:
		return new AdvectionDiffusionFactory(configuration.advection_diffusion_3D, _mesh);
	case PHYSICS::STRUCTURAL_MECHANICS_2D:
		return new StructuralMechanicsFactory(configuration.structural_mechanics_2D, _mesh);
	case PHYSICS::STRUCTURAL_MECHANICS_3D:
		return new StructuralMechanicsFactory(configuration.structural_mechanics_3D, _mesh);
	default:
		ESINFO(GLOBAL_ERROR) << "Unknown PHYSICS in configuration file";
		return NULL;
	}
}

template <class TType>
static void clear(std::vector<TType> &vector)
{
	for (size_t i = 0; i < vector.size(); i++) {
		delete vector[i];
	}
}

FactoryLoader::~FactoryLoader()
{
	clear(_instances);
	clear(_physics);
	clear(_linearSolvers);
	clear(_assemblers);
	clear(_timeStepSolvers);
	clear(_loadStepSolvers);
}

LinearSolver* FactoryLoader::getLinearSolver(const LoadStepSettingsBase &settings, Instance *instance) const
{
	switch (settings.solver_library) {
	case SOLVER_LIBRARY::ESPRESO:
		return new FETISolver(instance, settings.espreso);
	default:
		ESINFO(GLOBAL_ERROR) << "Not implemented requested SOLVER_LIBRARY.";
		return NULL;
	}
}

void Factory::loadPhysics(const GlobalConfiguration &configuration)
{
	_loader = createFactoryLoader(configuration);

	for (size_t step = 0; step < _loader->loadSteps(); step++) {
		_loadSteps.push_back(_loader->getLoadStepSolver(step, _mesh, _storeList));
	}

	_loader->preprocessMesh();
}

void FactoryLoader::preprocessMesh()
{
	// TODO: generalize it !!

	for (size_t i = 0; i < _physics.size(); i++) {

		switch (dynamic_cast<FETISolver*>(_linearSolvers.front())->configuration.method) {
		case ESPRESO_METHOD::TOTAL_FETI:
			_physics[i]->prepare();
			break;
		case ESPRESO_METHOD::HYBRID_FETI:
			switch (dynamic_cast<FETISolver*>(_linearSolvers.front())->configuration.B0_type) {
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
}

void Factory::finalize()
{
	// Detele store while finalizing because of Catalyst
	_storeList->finalize();
	_dispatcher.finalize();

	if (!_asyncStore) {
		// No store was setup, we need to delete the _asyncStore by ourself
		delete _asyncStore;
	}

}

Factory::~Factory()
{
	delete _mesh;
	delete _loader;
	delete _storeList;
}

void FactoryLoader::printError(const std::string &error) const
{
	ESINFO(GLOBAL_ERROR) << error;
}

}


