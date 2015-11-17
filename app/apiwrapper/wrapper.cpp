
#include "wrapper.h"

using namespace assembler;

std::list<FETI4IStructDoubleVector*> DataHolder::doubleVectors;
std::list<FETI4IStructIntVector*> DataHolder::intVectors;
std::list<FETI4IStructMap*> DataHolder::maps;
std::list<FETI4IStructMatrix*> DataHolder::matrices;
std::list<FETI4IStructFETIIntance*> DataHolder::instances;
MPI_Comm DataHolder::communicator;

using namespace assembler;

int FETI4IFinalize();

int FETI4IInit(MPI_Comm communicator)
{
	DataHolder::communicator = communicator;
	MPI_Comm_rank(communicator, &esconfig::MPIrank);
	MPI_Comm_size(communicator, &esconfig::MPIsize);
	return 0;
}

int FETI4ICreateMatrixElemental(
	FETI4IInt n,
	FETI4IInt nelt,
	FETI4IInt *eltptr,
	FETI4IInt *eltvar,
	double *values,
	FETI4IMatrix *stiffnessMatrix)
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
	FETI4IInt size,
	double *values,
	FETI4IDoubleVector *vector)
{
	DataHolder::doubleVectors.push_back(new FETI4IStructDoubleVector());
	DataHolder::doubleVectors.back()->size = size;
	DataHolder::doubleVectors.back()->values = new double[size];
	memcpy(DataHolder::doubleVectors.back()->values, values, size * sizeof(double));
	*vector = DataHolder::doubleVectors.back();
	return 0;
}

int FETI4ICreateIntVector(
	FETI4IInt size,
	FETI4IInt *values,
	FETI4IIntVector *vector)
{
	DataHolder::intVectors.push_back(new FETI4IStructIntVector());
	DataHolder::intVectors.back()->size = size;
	DataHolder::intVectors.back()->values = new FETI4IInt[size];
	memcpy(DataHolder::intVectors.back()->values, values, size * sizeof(FETI4IInt));
	*vector = DataHolder::intVectors.back();
	return 0;
}

int FETI4ICreateMap(
	FETI4IInt size,
	FETI4IInt *indices,
	double *values,
	FETI4IMap *vector)
{
	DataHolder::maps.push_back(new FETI4IStructMap());
	DataHolder::maps.back()->size = size;
	DataHolder::maps.back()->indices = new FETI4IInt[size];
	DataHolder::maps.back()->values = new double[size];
	memcpy(DataHolder::maps.back()->indices, indices, size * sizeof(FETI4IInt));
	memcpy(DataHolder::maps.back()->values, values, size * sizeof(double));
	*vector = DataHolder::maps.back();
	return 0;
}

int FETI4IPrepareFETIInstance(
	FETI4IMatrix *stiffnessMatrix,
	FETI4IDoubleVector *rhs,
	FETI4IMap *dirichlet,
	FETI4IIntVector *l2g,
	FETI4IIntVector *neighbourRanks,
	FETI4IFETIInstance *instance)
{
	API api((*stiffnessMatrix)->data, **rhs, **dirichlet, **l2g, **neighbourRanks);
	DataHolder::instances.push_back(new FETI4IStructFETIIntance(api));

	DataHolder::instances.back()->data.init();
	*instance = DataHolder::instances.back();
	return 0;
}

int FETI4ISolveFETI(
	FETI4IFETIInstance *instance,
	FETI4IInt size,
	double *values)
{
	std::vector<std::vector<double> > solution(1);
	solution[0] = std::vector<double>(values, values + size);
	(*instance)->data.solve(solution);
	memcpy(values, &solution[0][0], size);
	return 0;
}

int FETI4IDestroy(void *data)
{
	for (
		std::list<FETI4IStructDoubleVector*>::iterator it = DataHolder::doubleVectors.begin();
		it != DataHolder::doubleVectors.end();
		++it
	) {
		if (*it == data) {
			delete (*it)->values;
			delete *it;
			DataHolder::doubleVectors.erase(it);
			return 0;
		}
	}

	for (
		std::list<FETI4IStructIntVector*>::iterator it = DataHolder::intVectors.begin();
		it != DataHolder::intVectors.end();
		++it
	) {
		if (*it == data) {
			delete (*it)->values;
			delete *it;
			DataHolder::intVectors.erase(it);
			return 0;
		}
	}

	for (
		std::list<FETI4IStructMap*>::iterator it = DataHolder::maps.begin();
		it != DataHolder::maps.end();
		++it
	) {
		if (*it == data) {
			delete (*it)->indices;
			delete (*it)->values;
			delete *it;
			DataHolder::maps.erase(it);
			return 0;
		}
	}

	for (
		std::list<FETI4IStructMatrix*>::iterator it = DataHolder::matrices.begin();
		it != DataHolder::matrices.end();
		++it
	) {
		if (*it == data) {
			delete *it;
			DataHolder::matrices.erase(it);
			return 0;
		}
	}

	for (
		std::list<FETI4IStructFETIIntance*>::iterator it = DataHolder::instances.begin();
		it != DataHolder::instances.end();
		++it
	) {
		if (*it == data) {
			delete *it;
			DataHolder::instances.erase(it);
			return 0;
		}
	}

	return 0;
}

int FETI4IFinalize()
{
	return 0;
}


