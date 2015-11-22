
#include "wrapper.h"

using namespace assembler;

std::list<FETI4IStructRealVector*> DataHolder::doubleVectors;
std::list<FETI4IStructIntVector*> DataHolder::intVectors;
std::list<FETI4IStructMatrix*> DataHolder::matrices;
std::list<FETI4IStructIntance*> DataHolder::instances;

using namespace assembler;

int FETI4ICreateMatrixElemental(
		FETI4IMatrix *stiffnessMatrix,
		FETI4IInt n,
		FETI4IInt nelt,
		FETI4IInt* eltptr,
		FETI4IInt* eltvar,
		FETI4IReal* values)
{
	FETI4IInt indexing = eltptr[0];

	SparseVVPMatrix<eslocal> matrix(n, n);

	FETI4IInt value = 0;
	for (FETI4IInt e = 0; e < nelt; e++) {
		for (FETI4IInt i = eltptr[e] - indexing; i < eltptr[e + 1] - indexing; i++) {
			FETI4IInt size = eltptr[e + 1] - eltptr[e];
			for (FETI4IInt j = 0; j < size; j++) {
				matrix(eltvar[eltptr[e] - indexing + j] - indexing, eltvar[i] - indexing) = values[value++];
			}
		}
	}

	DataHolder::matrices.push_back(new FETI4IStructMatrix());
	DataHolder::matrices.back()->data = matrix;
	*stiffnessMatrix = DataHolder::matrices.back();
	return 0;
}

int FETI4ICreateDoubleVector(
		FETI4IRealVector *vector,
		FETI4IInt size,
		FETI4IReal* values)
{
	DataHolder::doubleVectors.push_back(new FETI4IStructRealVector());
	DataHolder::doubleVectors.back()->data = std::vector<double>(values, values + size);
	*vector = DataHolder::doubleVectors.back();
	return 0;
}

int FETI4ICreateIntVector(
		FETI4IIntVector *vector,
		FETI4IInt size,
		FETI4IInt *values)
{
	DataHolder::intVectors.push_back(new FETI4IStructIntVector());
	DataHolder::intVectors.back()->data = std::vector<eslocal>(values, values + size);
	*vector = DataHolder::intVectors.back();
	return 0;
}

int FETI4ICreateInstance(
		FETI4IInstance *instance,
		FETI4IInt* settings,
		FETI4IMatrix stiffnessMatrix,
		FETI4IRealVector rhs,
		FETI4IIntVector dirichlet_indices,
		FETI4IRealVector dirichlet_values,
		FETI4IIntVector l2g,
		FETI4IIntVector neighbourRanks,
		MPI_Comm communicator)
{
	MPI_Comm_rank(communicator, &esconfig::MPIrank);
	MPI_Comm_size(communicator, &esconfig::MPIsize);
	API api(
			stiffnessMatrix->data, rhs->data, dirichlet_indices->data, dirichlet_values->data,
			l2g->data, neighbourRanks->data);
	DataHolder::instances.push_back(new FETI4IStructIntance(api));

	DataHolder::instances.back()->data.init();
	*instance = DataHolder::instances.back();
	return 0;
}

int FETI4ISolve(
	FETI4IInstance instance,
	FETI4IInt solution_size,
	FETI4IReal* solution)
{
	std::vector<std::vector<double> > solutions(1);
	solutions[0] = std::vector<double>(solution, solution + solution_size);
	instance->data.solve(solutions);
	memcpy(solution, &solutions[0][0], solution_size);
	return 0;
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

int FETI4IDestroy(void *data)
{
	destroy(DataHolder::doubleVectors, data);
	destroy(DataHolder::intVectors, data);
	destroy(DataHolder::matrices, data);
	destroy(DataHolder::instances, data);

	return 0;
}


