
#include "wrapper.h"
#include "../../../libespreso/feti4i.h"

#include "../../assembler/physics/precomputed.h"
#include "../../assembler/physicssolver/assembler.h"
#include "../../assembler/physicssolver/timestep/linear.h"
#include "../../assembler/physicssolver/loadstep/steadystate.h"
#include "../../assembler/step.h"
#include "../../assembler/solution.h"

#include "../../output/resultstorelist.h"

#include "../../configuration/input/input.h"
#include "../../configuration/globalconfiguration.h"
#include "../../configuration/output.h"
#include "../../input/api/api.h"

#include "../../mesh/structures/mesh.h"
#include "../../solver/generic/FETISolver.h"

espreso::Environment espreso::DataHolder::environment;
std::list<FETI4IStructMatrix*> espreso::DataHolder::matrices;
std::list<FETI4IStructInstance*> espreso::DataHolder::instances;
espreso::TimeEval espreso::DataHolder::timeStatistics("API total time");

using namespace espreso;

FETI4IStructInstance::FETI4IStructInstance(FETI4IStructMatrix &matrix, eslocal *l2g, size_t size)
: instance(NULL), physics(NULL), linearSolver(NULL), assembler(NULL), timeStepSolver(NULL), loadStepSolver(NULL)
{
	output = new OutputConfiguration;
	store = new ResultStoreList(*output);
	mesh = new APIMesh(l2g, size);
	configuration = new ESPRESOSolver();
}

FETI4IStructInstance::~FETI4IStructInstance()
{
	if (instance != NULL) { delete instance; }
	if (physics != NULL) { delete physics; }
	if (linearSolver != NULL) { delete linearSolver; }
	if (assembler != NULL) { delete assembler; }
	if (timeStepSolver != NULL) { delete timeStepSolver; }
	if (loadStepSolver != NULL) { delete loadStepSolver; }

	delete store;
	delete output;
	delete mesh;
	delete configuration;
}

void FETI4ISetDefaultIntegerOptions(FETI4IInt* options)
{
	std::ifstream is("espreso.ecf");
	if (is.good()) {
		espreso::GlobalConfiguration configuration("espreso.ecf");

		options[FETI4I_SUBDOMAINS] = configuration.esdata.domains;

		options[FETI4I_ITERATIONS] = configuration.structural_mechanics_3D.physics_solver.load_steps_settings.at(1)->espreso.iterations;
		options[FETI4I_FETI_METHOD] = static_cast<int>(configuration.structural_mechanics_3D.physics_solver.load_steps_settings.at(1)->espreso.method);
		options[FETI4I_PRECONDITIONER] = static_cast<int>(configuration.structural_mechanics_3D.physics_solver.load_steps_settings.at(1)->espreso.preconditioner);
		options[FETI4I_CGSOLVER] = static_cast<int>(configuration.structural_mechanics_3D.physics_solver.load_steps_settings.at(1)->espreso.solver);
		options[FETI4I_N_MICS] = configuration.structural_mechanics_3D.physics_solver.load_steps_settings.at(1)->espreso.N_MICS;
		options[FETI4I_SC_SIZE] = configuration.structural_mechanics_3D.physics_solver.load_steps_settings.at(1)->espreso.SC_SIZE;

		options[FETI4I_VERBOSE_LEVEL] = environment->verbose_level;
		options[FETI4I_TESTING_LEVEL] = environment->testing_level;
		options[FETI4I_MEASURE_LEVEL] = environment->measure_level;
		options[FETI4I_PRINT_MATRICES] = environment->print_matrices;
		environment = &DataHolder::environment;
	} else {
		ESPRESOInput input;
		ESPRESOSolver solver;
		options[FETI4I_SUBDOMAINS] = input.domains;

		options[FETI4I_ITERATIONS] = solver.iterations;
		options[FETI4I_FETI_METHOD] = static_cast<int>(solver.method);
		options[FETI4I_PRECONDITIONER] = static_cast<int>(solver.preconditioner);
		options[FETI4I_CGSOLVER] = static_cast<int>(solver.solver);
		options[FETI4I_N_MICS] = solver.N_MICS;

		options[FETI4I_VERBOSE_LEVEL] = environment->verbose_level;
		options[FETI4I_TESTING_LEVEL] = environment->testing_level;
		options[FETI4I_MEASURE_LEVEL] = environment->measure_level;
		options[FETI4I_PRINT_MATRICES] = environment->print_matrices;

		options[FETI4I_SC_SIZE] = solver.SC_SIZE;
	}
}

void FETI4ISetDefaultRealOptions(FETI4IReal* options)
{
	std::ifstream is("espreso.ecf");
	if (is.good()) {
		espreso::GlobalConfiguration configuration("espreso.ecf");

		options[FETI4I_EPSILON] = configuration.structural_mechanics_3D.physics_solver.load_steps_settings.at(1)->espreso.epsilon;
		environment = &DataHolder::environment;
	} else {
		ESPRESOSolver solver;

		options[FETI4I_EPSILON] = solver.epsilon;
	}
}

static void FETI4ISetIntegerOptions(espreso::ESPRESOInput &input, espreso::ESPRESOSolver &solver, FETI4IInt* options)
{
	if (!input.parameters["domains"]->set(std::to_string(options[FETI4I_SUBDOMAINS]))) {
		ESINFO(GLOBAL_ERROR) << "Cannot set parameter 'SUBDOMAINS' to " << options[FETI4I_SUBDOMAINS];
	}
	if (!solver.parameters["iterations"]->set(std::to_string(options[FETI4I_ITERATIONS]))) {
		ESINFO(GLOBAL_ERROR) << "Cannot set parameter 'ITERATIONS' to " << options[FETI4I_ITERATIONS];
	}
	if (!solver.parameters["method"]->set(std::to_string(options[FETI4I_FETI_METHOD]))) {
		ESINFO(GLOBAL_ERROR) << "Cannot set parameter 'FETI_METHOD' to " << options[FETI4I_FETI_METHOD];
	}
	if (!solver.parameters["preconditioner"]->set(std::to_string(options[FETI4I_PRECONDITIONER]))) {
		ESINFO(GLOBAL_ERROR) << "Cannot set parameter 'PRECONDITIONER' to " << options[FETI4I_PRECONDITIONER];
	}
	if (!solver.parameters["solver"]->set(std::to_string(options[FETI4I_CGSOLVER]))) {
		ESINFO(GLOBAL_ERROR) << "Cannot set parameter 'CGSOLVER' to " << options[FETI4I_CGSOLVER];
	}
	if (!solver.parameters["N_MICS"]->set(std::to_string(options[FETI4I_N_MICS]))) {
		ESINFO(GLOBAL_ERROR) << "Cannot set parameter 'N_MICS' to " << options[FETI4I_N_MICS];
	}
	if (!solver.parameters["sc_size"]->set(std::to_string(options[FETI4I_SC_SIZE]))) {
		ESINFO(GLOBAL_ERROR) << "Cannot set parameter 'SC_SIZE' to " << options[FETI4I_SC_SIZE];
	}
	if (!environment->parameters["verbose_level"]->set(std::to_string(options[FETI4I_VERBOSE_LEVEL]))) {
		ESINFO(GLOBAL_ERROR) << "Cannot set parameter 'VERBOSE_LEVEL' to " << options[FETI4I_VERBOSE_LEVEL];
	}
	if (!environment->parameters["testing_level"]->set(std::to_string(options[FETI4I_TESTING_LEVEL]))) {
		ESINFO(GLOBAL_ERROR) << "Cannot set parameter 'TESTING_LEVEL' to " << options[FETI4I_TESTING_LEVEL];
	}
	if (!environment->parameters["measure_level"]->set(std::to_string(options[FETI4I_MEASURE_LEVEL]))) {
		ESINFO(GLOBAL_ERROR) << "Cannot set parameter 'MEASURE_LEVEL' to " << options[FETI4I_MEASURE_LEVEL];
	}
	if (!environment->parameters["print_matrices"]->set(std::to_string(options[FETI4I_PRINT_MATRICES]))) {
		ESINFO(GLOBAL_ERROR) << "Cannot set parameter 'PRINT_MATRICES' to " << options[FETI4I_PRINT_MATRICES];
	}

	solver.regularization = REGULARIZATION::NULL_PIVOTS;
	OutputConfiguration output;
	Reader::set(*environment, output);
}

static void FETI4ISetRealOptions(espreso::ESPRESOSolver &solver, FETI4IReal* options)
{
	solver.epsilon = options[FETI4I_EPSILON];
}

void FETI4ICreateStiffnessMatrix(
		FETI4IMatrix 	*matrix,
		FETI4IInt		type,
		FETI4IInt		indexBase)
{
	MPI_Comm_rank(environment->MPICommunicator, &environment->MPIrank);
	MPI_Comm_size(environment->MPICommunicator, &environment->MPIsize);
	OutputConfiguration output;
	Reader::set(*environment, output);

	DataHolder::timeStatistics.totalTime.startWithBarrier();
	TimeEvent event("Add element");
	DataHolder::timeStatistics.addEvent(event);

	DataHolder::matrices.push_back(new FETI4IStructMatrix(type, indexBase));
	*matrix = DataHolder::matrices.back();
}

void FETI4IAddElement(
		FETI4IMatrix	matrix,
		FETI4IInt		type,
		FETI4IInt		nodesSize,
		FETI4IInt*		nodes,
		FETI4IInt		dofsSize,
		FETI4IInt*		dofs,
		FETI4IReal*		values)
{
	espreso::DataHolder::timeStatistics.timeEvents.back().startWithoutBarrier();

	if (std::all_of(values, values + dofsSize, [] (double &value) { return value == 0; })) {
		// Skip elements with zero values
		return;
	}

	eslocal offset = matrix->offset;
	matrix->eType.push_back(type);
	matrix->eNodes.push_back(std::vector<eslocal>(nodes, nodes + nodesSize));
	matrix->eDOFs.push_back(std::vector<eslocal>(dofs, dofs + dofsSize));
	std::for_each(matrix->eNodes.back().begin(), matrix->eNodes.back().end(), [ &offset ] (eslocal &index) { index -= offset; });
	std::for_each(matrix->eDOFs.back().begin(), matrix->eDOFs.back().end(), [ &offset ] (eslocal &index) { index -= offset; });
	matrix->eMatrices.push_back(std::vector<double>(values, values + dofsSize * dofsSize));

	espreso::DataHolder::timeStatistics.timeEvents.back().endWithoutBarrier();
}

void FETI4ICreateInstance(
		FETI4IInstance 	*instance,
		FETI4IMatrix 	matrix,
		FETI4IInt 		size,
		FETI4IReal* 	rhs,
		FETI4IInt* 		l2g,
		FETI4IMPIInt 	neighbours_size,
		FETI4IMPIInt*	neighbours,
		FETI4IInt 		dirichlet_size,
		FETI4IInt* 		dirichlet_indices,
		FETI4IReal* 	dirichlet_values,
		FETI4IInt* 		integer_options,
		FETI4IReal*		real_options)
{
	DataHolder::instances.push_back(new FETI4IStructInstance(*matrix, l2g, size));

	ESPRESOInput input;
	FETI4ISetIntegerOptions(input, *DataHolder::instances.back()->configuration, integer_options);
	FETI4ISetRealOptions(*DataHolder::instances.back()->configuration, real_options);

	TimeEvent event("Create FETI4I instance"); event.startWithBarrier();

	ESINFO(OVERVIEW) << "ESPRESO create solver instance";

	std::vector<int> neighClusters = std::vector<int>(neighbours, neighbours + neighbours_size);

	input::API::load(
			input,
			*DataHolder::instances.back()->mesh, matrix->offset,
			matrix->eType, matrix->eNodes, matrix->eDOFs, matrix->eMatrices,
			dirichlet_size, dirichlet_indices, dirichlet_values,
			neighClusters,
			size, l2g);

	*instance = DataHolder::instances.back();

	DataHolder::instances.back()->instance = new Instance(*DataHolder::instances.back()->mesh);
	DataHolder::instances.back()->physics = new Precomputed(DataHolder::instances.back()->mesh, DataHolder::instances.back()->instance, (espreso::MatrixType)matrix->type, rhs, size);
	DataHolder::instances.back()->linearSolver = new FETISolver(DataHolder::instances.back()->instance, *DataHolder::instances.back()->configuration);
	DataHolder::instances.back()->assembler = new Assembler(
			*DataHolder::instances.back()->instance,
			*DataHolder::instances.back()->physics,
			*DataHolder::instances.back()->mesh,
			*DataHolder::instances.back()->store,
			*DataHolder::instances.back()->linearSolver);
	DataHolder::instances.back()->timeStepSolver = new LinearTimeStep(*DataHolder::instances.back()->assembler);
	DataHolder::instances.back()->loadStepSolver = new SteadyStateSolver(*DataHolder::instances.back()->timeStepSolver, 1);

	switch (DataHolder::instances.back()->configuration->method) {
	case ESPRESO_METHOD::TOTAL_FETI:
		DataHolder::instances.back()->physics->prepare();
		break;
	case ESPRESO_METHOD::HYBRID_FETI:
		DataHolder::instances.back()->physics->prepareHybridTotalFETIWithKernels();
		break;
	default:
		ESINFO(ERROR) << "API request unknown FETI method.";
	}

	event.endWithBarrier(); DataHolder::timeStatistics.addEvent(event);
}

void FETI4ISolve(
		FETI4IInstance 	instance,
		FETI4IInt 		solution_size,
		FETI4IReal*		solution)
{
	TimeEvent event("Solve FETI4I instance"); event.startWithBarrier();

	Step step;
	Logging::step = &step;
	instance->loadStepSolver->run(step);

	memcpy(solution, instance->instance->solutions[espreso::Precomputed::SolutionIndex::MERGED]->data[0].data(), solution_size * sizeof(double));

	event.endWithBarrier(); DataHolder::timeStatistics.addEvent(event);
	DataHolder::timeStatistics.totalTime.endWithBarrier();
	DataHolder::timeStatistics.printStatsMPI();
}

template <typename TFETI4I>
static void destroy(std::list<TFETI4I*> &list, void *value)
{
	for (typename std::list<TFETI4I*>::iterator it = list.begin(); it != list.end(); ++it) {
		if (*it == value) {
			delete *it;
			list.erase(it);
			return;
		}
	}
}

void FETI4IDestroy(void *data)
{
	destroy(DataHolder::matrices, data);
	destroy(DataHolder::instances, data);
}

