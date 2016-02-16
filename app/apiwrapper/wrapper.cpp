
#include "wrapper.h"

using namespace assembler;
//using namespace esinput;

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
		ESLOG(eslog::ERROR) << "Cannot read file " << fileName;
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
		ESLOG(eslog::ERROR) << "Cannot read file " << fileName;
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
//	MPI_Comm_rank(MPI_COMM_WORLD, &esconfig::MPIrank);
//	MPI_Comm_size(MPI_COMM_WORLD, &esconfig::MPIsize);
//
//	CubeSettings cube(esconfig::MPIrank, esconfig::MPIsize);
//	MeshGenerator generator(new CubeGenerator<Hexahedron8>(cube));
//
//	mesh::Mesh mesh(esconfig::MPIrank, esconfig::MPIsize);
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

void FETI4ICreateStiffnessMatrix(
		FETI4IMatrix 	*matrix,
		FETI4IInt		indexBase)
{
	DataHolder::matrices.push_back(new FETI4IStructMatrix(indexBase));
	DataHolder::matrices.back()->K.resize(esconfig::mesh::subdomains, SparseVVPMatrix<eslocal>(0, 0));
	*matrix = DataHolder::matrices.back();
}

void FETI4IAddElement(
		FETI4IMatrix 	matrix,
		FETI4IInt 		size,
		FETI4IInt* 		indices,
		FETI4IReal* 	values)
{
	eslocal offset = matrix->offset;
	if (esconfig::mesh::subdomains == 1) {
		for (eslocal i = 0; i < size; i++) {
			for (eslocal j = 0; j < size; j++) {
				matrix->K[0](indices[i] - offset,  indices[j] - offset) = values[i * size + j];
			}
		}
	} else {
		matrix->eIndices.push_back(std::vector<eslocal>(indices, indices + size));
		std::for_each(matrix->eIndices.back().begin(), matrix->eIndices.back().end(), [ &offset ] (eslocal &index) { index -= offset; });
		matrix->eMatrix.push_back(std::vector<double>(values, values + size));
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

	API api;
	DataHolder::instances.push_back(new FETI4IStructInstance(api));
	DataHolder::instances.back()->K.resize(matrix->K.size(), SparseCSRMatrix<eslocal>(0, 0));
	for (size_t i = 0; i < matrix->K.size(); i++) {
		DataHolder::instances.back()->K[i] = matrix->K[i];
	}
	api.K = &(DataHolder::instances.back()->K);
	api.indexing = matrix->offset;
	api.size = size;
	api.rhs = rhs;
	for (size_t i = 0; i < dirichlet_size; i++) {
		dirichlet_indices[i] = l2g[dirichlet_indices[i] - matrix->offset];
	}

	api.dirichlet_size = dirichlet_size;
	api.dirichlet_indices = dirichlet_indices;
	api.dirichlet_values = dirichlet_values;
	api.l2g = l2g;
	api.neighbours_size = neighbours_size;
	api.neighbours = neighbours;
	DataHolder::instances.back()->data = assembler::LinearElasticity<assembler::API>(api);

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
	memcpy(solution, &solutions[0][0], solution_size * sizeof(double));

	instance->data.finalize();
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
	MPI_Comm_rank(MPI_COMM_WORLD, &esconfig::MPIrank);
	MPI_Comm_size(MPI_COMM_WORLD, &esconfig::MPIsize);
	std::vector<FETI4IInt> eCount;
	std::vector<FETI4IInt> eSize;
	std::stringstream ssEl, ssKi;
	ssEl << "examples/api/cube/" << esconfig::MPIrank << "/elements.txt";
	ssKi << "examples/api/cube/" << esconfig::MPIrank << "/Ki0.txt";
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
	ssKi << "examples/api/cube/" << esconfig::MPIrank << "/Ki" << index << ".txt";
	ssKv << "examples/api/cube/" << esconfig::MPIrank << "/Ke" << index << ".bin";
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
	ssRhs << "examples/api/cube/" << esconfig::MPIrank << "/rhs.txt";
	ssD << "examples/api/cube/" << esconfig::MPIrank << "/dirichlet_indices.txt";
	ssN << "examples/api/cube/" << esconfig::MPIrank << "/neighbours.txt";
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
	ssRhs << "examples/api/cube/" << esconfig::MPIrank << "/rhs.txt";
	ssL << "examples/api/cube/" << esconfig::MPIrank << "/l2g.txt";
	ssDi << "examples/api/cube/" << esconfig::MPIrank << "/dirichlet_indices.txt";
	ssDv << "examples/api/cube/" << esconfig::MPIrank << "/dirichlet_values.txt";
	ssN << "examples/api/cube/" << esconfig::MPIrank << "/neighbours.txt";
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

