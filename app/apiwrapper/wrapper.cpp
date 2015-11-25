
#include "wrapper.h"

using namespace assembler;
using namespace esinput;

std::list<FETI4IStructMatrix*> DataHolder::matrices;
std::list<FETI4IStructInstance*> DataHolder::instances;

using namespace assembler;

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


