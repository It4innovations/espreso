
#include "wrapper.h"

using namespace assembler;
using namespace esinput;

std::list<FETI4IStructMatrix*> DataHolder::matrices;
std::list<FETI4IStructInstance*> DataHolder::instances;

using namespace assembler;

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
		std::cerr << "Cannot read file " << fileName << "\n";
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
		std::cerr << "Cannot read file " << fileName << "\n";
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
	MPI_Comm_rank(MPI_COMM_WORLD, &esconfig::MPIrank);
	MPI_Comm_size(MPI_COMM_WORLD, &esconfig::MPIsize);

	CubeSettings cube(esconfig::MPIrank, esconfig::MPIsize);
	MeshGenerator generator(new CubeGenerator<Hexahedron8>(cube));

	mesh::Mesh mesh(esconfig::MPIrank, esconfig::MPIsize);
	generator.load(mesh);

	FEM fem(mesh);
	LinearElasticity<FEM> solver(fem);

	std::vector<std::vector<double> > solution;

	solver.init();
	solver.solve(solution);
}

void FETI4ICreateStiffnessMatrix(
		FETI4IMatrix 	*matrix,
		FETI4IInt		indexBase)
{
	DataHolder::matrices.push_back(new FETI4IStructMatrix(indexBase));
	*matrix = DataHolder::matrices.back();
}

void FETI4IAddElement(
		FETI4IMatrix 	matrix,
		FETI4IInt 		size,
		FETI4IInt* 		indices,
		FETI4IReal* 	values)
{
	eslocal offset = matrix->offset;
	for (eslocal i = 0; i < size; i++) {
		for (eslocal j = 0; j < size; j++) {
			matrix->data(indices[i] - offset,  indices[j] - offset) = values[i * size + j];
		}
	}
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
		FETI4IInt* 		dirichlet_indices,  //TODO which numbering? we prefer global numbering
		FETI4IReal* 	dirichlet_values)
{
	MPI_Comm_rank(MPI_COMM_WORLD, &esconfig::MPIrank);
	MPI_Comm_size(MPI_COMM_WORLD, &esconfig::MPIsize);

	API2 api;
	DataHolder::instances.push_back(new FETI4IStructInstance(api));
	DataHolder::instances.back()->K = matrix->data;
	api.K = &(DataHolder::instances.back()->K);
	api.size = size;
	api.rhs = rhs;
	api.dirichlet_size = dirichlet_size;
	api.dirichlet_indices = dirichlet_indices;
	api.dirichlet_values = dirichlet_values;
	api.l2g = l2g;
	api.neighbours_size = neighbours_size;
	api.neighbours = neighbours;
	DataHolder::instances.back()->data = assembler::LinearElasticity<assembler::API2>(api);

	DataHolder::instances.back()->data.init();
	*instance = DataHolder::instances.back();
}

void FETI4ISolve(
		FETI4IInstance 	instance,
		FETI4IInt 		solution_size,
		FETI4IReal*		solution)
{
	std::vector<std::vector<double> > solutions(1);
	solutions[0] = std::vector<double>(solution, solution + solution_size);
	instance->data.solve(solutions);
	memcpy(solution, &solutions[0][0], solution_size);
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
	std::vector<FETI4IInt> eCount;
	std::vector<FETI4IInt> eSize;
	readFile(eCount, "examples/api/cube/0/elements.txt");
	readFile(eSize, "examples/api/cube/0/Ki.txt");

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
	ssKi << "examples/api/cube/0/Ki" << index << ".txt";
	ssKv << "examples/api/cube/0/Ke" << index << ".bin";
	readFile(Ki, ssKi.str().c_str());
	readBinary(Kv, ssKv.str().c_str());
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

	readFile(rhs, "examples/api/cube/0/rhs.txt");
	readFile(dirichlet, "examples/api/cube/0/dirichlet_indices.txt");
	readFile(neighbours, "examples/api/cube/0/neighbours.txt");

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
	readFile(_rhs, "examples/api/cube/0/rhs.txt");
	readFile(_dirichlet_indices, "examples/api/cube/0/dirichlet_indices.txt");
	readFile(_dirichlet_values, "examples/api/cube/0/dirichlet_values.txt");
	readFile(_l2g, "examples/api/cube/0/l2g.txt");
	readFile(_neighbours, "examples/api/cube/0/neighbours.txt");

	memcpy(rhs, _rhs.data(), _rhs.size() * sizeof(double));
	memcpy(l2g, _l2g.data(), _l2g.size() * sizeof(eslocal));
	memcpy(dirichlet_indices, _dirichlet_indices.data(), _dirichlet_indices.size() * sizeof(eslocal));
	memcpy(dirichlet_values, _dirichlet_values.data(), _dirichlet_values.size() * sizeof(double));
	memcpy(neighbours, _neighbours.data(), _neighbours.size() * sizeof(eslocal));
}

