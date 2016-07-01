
#include "wrapper.h"

std::list<FETI4IStructMatrix*> espreso::DataHolder::matrices;
std::list<FETI4IStructInstance*> espreso::DataHolder::instances;
espreso::TimeEval espreso::DataHolder::timeStatistics("API total time");

using namespace espreso;

template <typename Ttype>
static void readFile(typename std::vector<Ttype> &vector, std::string fileName) {
	vector.clear();
	std::ifstream file(fileName.c_str());
	if (file.is_open()) {
		Ttype value;
		while (file >> value) {
			vector.push_back(value);
		}
	} else {
		ESINFO(ERROR) << "Cannot read file " << fileName;
		exit(EXIT_FAILURE);
	}
}

static void readBinary(std::vector<double> &vector, std::string fileName) {
	std::ifstream file(fileName.c_str(), std::fstream::binary);
	if (file.is_open()) {
		for (size_t i = 0; i < vector.size(); i++) {
			double value;
			file.read(reinterpret_cast<char *>(&value), sizeof(double));
			vector[i] = value;
		}
	} else {
		ESINFO(ERROR) << "Cannot read file " << fileName;
		exit(EXIT_FAILURE);
	}
}

//int FETI4ICreateStiffnessMatrix(
//		FETI4IMatrix *stiffnessMatrix,
//		FETI4IInt n,
//		FETI4IInt nelt,
//		FETI4IInt* eltptr,
//		FETI4IInt* eltvar,
//		FETI4IReal* values)
//{
//	FETI4IInt indexing = eltptr[0];
//
//	SparseVVPMatrix<eslocal> matrix(n, n);
//
//	FETI4IInt value = 0;
//	for (FETI4IInt e = 0; e < nelt; e++) {
//		for (FETI4IInt i = eltptr[e] - indexing; i < eltptr[e + 1] - indexing; i++) {
//			FETI4IInt size = eltptr[e + 1] - eltptr[e];
//			for (FETI4IInt j = 0; j < size; j++) {
//				matrix(eltvar[eltptr[e] - indexing + j] - indexing, eltvar[i] - indexing) = values[value++];
//			}
//		}
//	}
//
//	DataHolder::matrices.push_back(new FETI4IStructMatrix(indexing));
//	DataHolder::matrices.back()->data = matrix;
//	*stiffnessMatrix = DataHolder::matrices.back();
//	return 0;
//}

void FETI4ITest()
{
//	MPI_Comm_rank(MPI_COMM_WORLD, &esconfig::env::MPIrank);
//	MPI_Comm_size(MPI_COMM_WORLD, &esconfig::env::MPIsize);
//
//	CubeSettings cube(esconfig::env::MPIrank, esconfig::env::MPIsize);
//	MeshGenerator generator(new CubeGenerator<Hexahedron8>(cube));
//
//	mesh::Mesh mesh(esconfig::env::MPIrank, esconfig::env::MPIsize);
//	generator.load(mesh);
//
//	FEM fem(mesh);
//	LinearElasticity<FEM> solver(fem);
//
//	std::vector<std::vector<double> > solution;
//
//	solver.init();
//	solver.solve(solution);
}


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
		FETI4IMatrix 	matrix,
		FETI4IInt 		size,
		FETI4IInt* 		indices,
		FETI4IReal* 	values)
{
	espreso::DataHolder::timeStatistics.timeEvents.back().startWithoutBarrier();

	if (std::all_of(values, values + size, [] (double &value) { return value == 0; })) {
		// Skip elements with zero values
		return;
	}

	eslocal offset = matrix->offset;
	matrix->eIndices.push_back(std::vector<eslocal>(indices, indices + size));
	std::for_each(matrix->eIndices.back().begin(), matrix->eIndices.back().end(), [ &offset ] (eslocal &index) { index -= offset; });
	matrix->eMatrices.push_back(std::vector<double>(values, values + size * size));

	espreso::DataHolder::timeStatistics.timeEvents.back().endWithoutBarrier();
}

void FETI4ICreateInstance(
		FETI4IInstance 	*instance,
		FETI4IMatrix 	matrix,
		FETI4IInt 		size,
		FETI4IReal* 	rhs,
		FETI4IInt* 		l2g,                     /* length of both rhs and l2g is size */
		FETI4IMPIInt 	neighbours_size,
		FETI4IMPIInt*	neighbours,
		FETI4IInt 		dirichlet_size,
		FETI4IInt* 		dirichlet_indices,
		FETI4IReal* 	dirichlet_values,
		FETI4IInt* 		integer_options,
		FETI4IReal* 		real_options)
{
	FETI4ISetIntegerOptions(integer_options);
	FETI4ISetRealOptions(real_options);

	TimeEvent event("Create FETI4I instance"); event.startWithBarrier();

	ESINFO(OVERVIEW) << "ESPRESO create solver instance";

	std::vector<eslocal> neighClusters = std::vector<eslocal>(neighbours, neighbours + neighbours_size);

	APIMesh *mesh = new APIMesh(matrix->eMatrices);
	input::API::load(*mesh, matrix->eIndices, neighClusters, size, l2g);

	API api(mesh);
	api.indexing = matrix->offset;
	api.size = size;
	api.rhs = rhs;
	api.dirichlet_size = dirichlet_size;
	api.dirichlet_indices = dirichlet_indices;
	api.dirichlet_values = dirichlet_values;
	api.l2g = l2g;
	api.neighbours_size = neighbours_size;
	api.neighbours = neighbours;

	DataHolder::instances.push_back(new FETI4IStructInstance(api, mesh));
	DataHolder::instances.back()->data.init();
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
	instance->data.solve(solutions);
	memcpy(solution, &solutions[0][0], solution_size * sizeof(double));

	instance->data.finalize();

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


void TEST4IGetElementsInfo(
		FETI4IInt		*elements,
		FETI4IInt		*elementSize)
{
	MPI_Comm_rank(MPI_COMM_WORLD, &config::env::MPIrank);
	MPI_Comm_size(MPI_COMM_WORLD, &config::env::MPIsize);
	std::vector<FETI4IInt> eCount;
	std::vector<FETI4IInt> eSize;
	std::stringstream ssEl, ssKi;
	ssEl << "examples/api/cube/" << config::env::MPIrank << "/elements.txt";
	ssKi << "examples/api/cube/" << config::env::MPIrank << "/Ki0.txt";
	readFile(eCount, ssEl.str());
	readFile(eSize, ssKi.str());

	*elements = eCount.size();
	*elementSize = eSize.size();
}

void TEST4IGetElement(
		FETI4IInt		index,
		FETI4IInt*		*indices,
		FETI4IReal*		*values)
{
	std::vector<FETI4IInt> Ki;
	std::vector<FETI4IReal> Kv;

	std::stringstream ssKi, ssKv;
	ssKi << "examples/api/cube/" << config::env::MPIrank << "/Ki" << index << ".txt";
	ssKv << "examples/api/cube/" << config::env::MPIrank << "/Ke" << index << ".bin";
	readFile(Ki, ssKi.str());
	Kv.resize(Ki.size() * Ki.size());
	readBinary(Kv, ssKv.str());
	memcpy(indices, Ki.data(), Ki.size() * sizeof(eslocal));
	memcpy(values, Kv.data(), Kv.size() * sizeof(double));
}

void TEST4IGetInstanceInfo(
		FETI4IInt		*rhs_size,
		FETI4IInt		*dirichlet_size,
		FETI4IInt		*neighbours_size)
{
	std::vector<FETI4IReal> rhs;
	std::vector<FETI4IInt> dirichlet;
	std::vector<FETI4IMPIInt> neighbours;

	std::stringstream ssRhs, ssD, ssN;
	ssRhs << "examples/api/cube/" << config::env::MPIrank << "/rhs.txt";
	ssD << "examples/api/cube/" << config::env::MPIrank << "/dirichlet_indices.txt";
	ssN << "examples/api/cube/" << config::env::MPIrank << "/neighbours.txt";
	readFile(rhs, ssRhs.str().c_str());
	readFile(dirichlet, ssD.str().c_str());
	readFile(neighbours, ssN.str().c_str());

	*rhs_size = rhs.size();
	*dirichlet_size = dirichlet.size();
	*neighbours_size = neighbours.size();
}

void TEST4IGetInstance(
		FETI4IReal*		*rhs,
		FETI4IInt*		*l2g,
		FETI4IInt*		*dirichlet_indices,
		FETI4IReal*		*dirichlet_values,
		FETI4IInt*		*neighbours)
{
	std::vector<FETI4IReal> _rhs;
	std::vector<FETI4IInt> _l2g;
	std::vector<FETI4IInt> _dirichlet_indices;
	std::vector<FETI4IReal> _dirichlet_values;
	std::vector<FETI4IMPIInt> _neighbours;

	std::stringstream ssRhs, ssL, ssDi, ssDv, ssN;
	ssRhs << "examples/api/cube/" << config::env::MPIrank << "/rhs.txt";
	ssL << "examples/api/cube/" << config::env::MPIrank << "/l2g.txt";
	ssDi << "examples/api/cube/" << config::env::MPIrank << "/dirichlet_indices.txt";
	ssDv << "examples/api/cube/" << config::env::MPIrank << "/dirichlet_values.txt";
	ssN << "examples/api/cube/" << config::env::MPIrank << "/neighbours.txt";
	readFile(_rhs, ssRhs.str());
	readFile(_dirichlet_indices, ssDi.str());
	readFile(_dirichlet_values, ssDv.str());
	readFile(_l2g, ssL.str());
	readFile(_neighbours, ssN.str());

	memcpy(rhs, _rhs.data(), _rhs.size() * sizeof(double));
	memcpy(l2g, _l2g.data(), _l2g.size() * sizeof(eslocal));
	memcpy(dirichlet_indices, _dirichlet_indices.data(), _dirichlet_indices.size() * sizeof(eslocal));
	memcpy(dirichlet_values, _dirichlet_values.data(), _dirichlet_values.size() * sizeof(double));
	memcpy(neighbours, _neighbours.data(), _neighbours.size() * sizeof(eslocal));
}
