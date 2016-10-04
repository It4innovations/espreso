
#include "wrapper.h"
#include "../../../libespreso/feti4i.h"
#include "../../assembler/instance/precomputed/instance.h"
#include "../../assembler/physics/precomputed/uniformsymmetric3dofs/assembler.h"

std::list<FETI4IStructMatrix*> espreso::DataHolder::matrices;
std::list<FETI4IStructInstance*> espreso::DataHolder::instances;
espreso::TimeEval espreso::DataHolder::timeStatistics("API total time");

using namespace espreso;

void FETI4ISetDefaultIntegerOptions(FETI4IInt* options)
{
	options[FETI4I_SUBDOMAINS] = config::mesh::SUBDOMAINS;

	options[FETI4I_ITERATIONS] = config::solver::ITERATIONS;
	options[FETI4I_FETI_METHOD] = static_cast<int>(config::solver::FETI_METHOD);
	options[FETI4I_PRECONDITIONER] = static_cast<int>(config::solver::PRECONDITIONER);
	options[FETI4I_CGSOLVER] = static_cast<int>(config::solver::CGSOLVER);
	options[FETI4I_N_MICS] = config::solver::N_MICS;

	options[FETI4I_VERBOSE_LEVEL] = config::info::VERBOSE_LEVEL;
	options[FETI4I_TESTING_LEVEL] = config::info::TESTING_LEVEL;
	options[FETI4I_MEASURE_LEVEL] = config::info::MEASURE_LEVEL;
	options[FETI4I_PRINT_MATRICES] = config::info::PRINT_MATRICES;
}

void FETI4ISetDefaultRealOptions(FETI4IReal* options)
{
	options[FETI4I_EPSILON] = config::solver::EPSILON;
}

static void FETI4ISetIntegerOptions(FETI4IInt* options)
{
	ParametersReader reader(config::parameters);

	if (!reader.setParameter(&config::mesh::SUBDOMAINS, options[FETI4I_SUBDOMAINS])) {
		ESINFO(GLOBAL_ERROR) << "Cannot set parameter 'SUBDOMAINS' to " << options[FETI4I_SUBDOMAINS];
	}
	if (!reader.setParameter(&config::solver::ITERATIONS, options[FETI4I_ITERATIONS])) {
		ESINFO(GLOBAL_ERROR) << "Cannot set parameter 'ITERATIONS' to " << options[FETI4I_ITERATIONS];
	}
	if (!reader.setParameter(&config::solver::FETI_METHOD, options[FETI4I_FETI_METHOD])) {
		ESINFO(GLOBAL_ERROR) << "Cannot set parameter 'FETI_METHOD' to " << options[FETI4I_FETI_METHOD];
	}
	if (!reader.setParameter(&config::solver::PRECONDITIONER, options[FETI4I_PRECONDITIONER])) {
		ESINFO(GLOBAL_ERROR) << "Cannot set parameter 'PRECONDITIONER' to " << options[FETI4I_PRECONDITIONER];
	}
	if (!reader.setParameter(&config::solver::CGSOLVER, options[FETI4I_CGSOLVER])) {
		ESINFO(GLOBAL_ERROR) << "Cannot set parameter 'CGSOLVER' to " << options[FETI4I_CGSOLVER];
	}
	if (!reader.setParameter(&config::solver::N_MICS, options[FETI4I_N_MICS])) {
		ESINFO(GLOBAL_ERROR) << "Cannot set parameter 'N_MICS' to " << options[FETI4I_N_MICS];
	}
	if (!reader.setParameter(&config::info::VERBOSE_LEVEL, options[FETI4I_VERBOSE_LEVEL])) {
		ESINFO(GLOBAL_ERROR) << "Cannot set parameter 'VERBOSE_LEVEL' to " << options[FETI4I_VERBOSE_LEVEL];
	}
	if (!reader.setParameter(&config::info::TESTING_LEVEL, options[FETI4I_TESTING_LEVEL])) {
		ESINFO(GLOBAL_ERROR) << "Cannot set parameter 'TESTING_LEVEL' to " << options[FETI4I_TESTING_LEVEL];
	}
	if (!reader.setParameter(&config::info::MEASURE_LEVEL, options[FETI4I_MEASURE_LEVEL])) {
		ESINFO(GLOBAL_ERROR) << "Cannot set parameter 'MEASURE_LEVEL' to " << options[FETI4I_MEASURE_LEVEL];
	}
	if (!reader.setParameter(&config::info::PRINT_MATRICES, options[FETI4I_PRINT_MATRICES])) {
		ESINFO(GLOBAL_ERROR) << "Cannot set parameter 'PRINT_MATRICES' to " << options[FETI4I_PRINT_MATRICES];
	}
}

static void FETI4ISetRealOptions(FETI4IReal* options)
{
	ParametersReader reader(config::parameters);

	if (!reader.setParameter(&config::solver::EPSILON, options[FETI4I_EPSILON])) {
		ESINFO(GLOBAL_ERROR) << "Cannot set parameter 'EPSILON' to " << options[FETI4I_PRINT_MATRICES];
	}
}

void FETI4ICreateStiffnessMatrix(
		FETI4IMatrix 	*matrix,
		FETI4IInt		indexBase)
{
	MPI_Comm_rank(MPI_COMM_WORLD, &config::env::MPIrank);
	MPI_Comm_size(MPI_COMM_WORLD, &config::env::MPIsize);
	config::solver::REGULARIZATION = config::solver::REGULARIZATIONalternative::NULL_PIVOTS;

	DataHolder::timeStatistics.totalTime.startWithBarrier();
	TimeEvent event("Add element");
	DataHolder::timeStatistics.addEvent(event);

	DataHolder::matrices.push_back(new FETI4IStructMatrix(indexBase));
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
	FETI4ISetIntegerOptions(integer_options);
	FETI4ISetRealOptions(real_options);
	ParametersReader::printParameters(config::parameters, config::info::VERBOSE_LEVEL);

	TimeEvent event("Create FETI4I instance"); event.startWithBarrier();

	ESINFO(OVERVIEW) << "ESPRESO create solver instance";

	std::vector<eslocal> neighClusters = std::vector<eslocal>(neighbours, neighbours + neighbours_size);

	DataHolder::instances.push_back(new FETI4IStructInstance(*matrix));
	input::API::load(
			DataHolder::instances.back()->mesh, matrix->offset,
			matrix->eType, matrix->eNodes, matrix->eDOFs,
			dirichlet_size, dirichlet_indices, dirichlet_values,
			neighClusters,
			size, l2g);

	DataHolder::instances.back()->instance = new PrecomputedInstance<EqualityConstraints, UniformSymmetric3DOFs>(DataHolder::instances.back()->mesh, rhs, size);
	DataHolder::instances.back()->instance->init();
	*instance = DataHolder::instances.back();

	event.endWithBarrier(); DataHolder::timeStatistics.addEvent(event);
}

void FETI4ISolve(
		FETI4IInstance 	instance,
		FETI4IInt 		solution_size,
		FETI4IReal*		solution)
{
	TimeEvent event("Solve FETI4I instance"); event.startWithBarrier();

	std::vector<std::vector<double> > solutions(1);
	solutions[0] = std::vector<double>(solution, solution + solution_size);
	instance->instance->solve(solutions);
	memcpy(solution, &solutions[0][0], solution_size * sizeof(double));

	instance->instance->finalize();

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

