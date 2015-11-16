
#include "wrapper.h"

using namespace assembler;

std::list<ESPRESOStructDoubleVector*> DataHolder::doubleVectors;
std::list<ESPRESOStructIntVector*> DataHolder::intVectors;
std::list<ESPRESOStructMap*> DataHolder::maps;
std::list<ESPRESOStructMatrix*> DataHolder::matrices;
std::list<ESPRESOStructFETIIntance*> DataHolder::instances;
MPI_Comm DataHolder::communicator;

using namespace assembler;

int ESPRESOFinalize();

int ESPRESOInit(MPI_Comm communicator)
{
	DataHolder::communicator = communicator;
	MPI_Comm_rank(communicator, &esconfig::MPIrank);
	MPI_Comm_size(communicator, &esconfig::MPIsize);
	return 0;
}

int ESPRESOCreateMatrixElemental(
	esint n,
	esint nelt,
	esint *eltptr,
	esint *eltvar,
	double *values,
	ESPRESOMatrix *stiffnessMatrix)
{
	esint indexing = eltptr[0];

	SparseVVPMatrix<eslocal> matrix(n, n);

	esint value = 0;
	for (esint e = 0; e < nelt; e++) {
		for (esint i = eltptr[e] - indexing; i < eltptr[e + 1] - indexing; i++) {
			esint size = eltptr[e + 1] - eltptr[e];
			for (esint j = 0; j < size; j++) {
				matrix(eltvar[eltptr[e] - indexing + j] - indexing, eltvar[i] - indexing) = values[value++];
			}
		}
	}

	DataHolder::matrices.push_back(new ESPRESOStructMatrix());
	DataHolder::matrices.back()->data = matrix;
	*stiffnessMatrix = DataHolder::matrices.back();
	return 0;
}

int ESPRESOCreateDoubleVector(
	esint size,
	double *values,
	ESPRESODoubleVector *vector)
{
	DataHolder::doubleVectors.push_back(new ESPRESOStructDoubleVector());
	DataHolder::doubleVectors.back()->size = size;
	DataHolder::doubleVectors.back()->values = new double[size];
	memcpy(DataHolder::doubleVectors.back()->values, values, size * sizeof(double));
	*vector = DataHolder::doubleVectors.back();
	return 0;
}

int ESPRESOCreateIntVector(
	esint size,
	esint *values,
	ESPRESOIntVector *vector)
{
	DataHolder::intVectors.push_back(new ESPRESOStructIntVector());
	DataHolder::intVectors.back()->size = size;
	DataHolder::intVectors.back()->values = new esint[size];
	memcpy(DataHolder::intVectors.back()->values, values, size * sizeof(esint));
	*vector = DataHolder::intVectors.back();
	return 0;
}

int ESPRESOCreateMap(
	esint size,
	esint *indices,
	double *values,
	ESPRESOMap *vector)
{
	DataHolder::maps.push_back(new ESPRESOStructMap());
	DataHolder::maps.back()->size = size;
	DataHolder::maps.back()->indices = new esint[size];
	DataHolder::maps.back()->values = new double[size];
	memcpy(DataHolder::maps.back()->indices, indices, size * sizeof(esint));
	memcpy(DataHolder::maps.back()->values, values, size * sizeof(double));
	*vector = DataHolder::maps.back();
	return 0;
}

int ESPRESOPrepareFETIInstance(
	ESPRESOMatrix *stiffnessMatrix,
	ESPRESODoubleVector *rhs,
	ESPRESOMap *dirichlet,
	ESPRESOIntVector *l2g,
	ESPRESOIntVector *neighbourRanks,
	ESPRESOFETIInstance *instance)
{
	API api((*stiffnessMatrix)->data, **rhs, **dirichlet, **l2g, **neighbourRanks);
	DataHolder::instances.push_back(new ESPRESOStructFETIIntance(api));

	DataHolder::instances.back()->data.init();
	*instance = DataHolder::instances.back();
	return 0;
}

int ESPRESOSolveFETI(
	ESPRESOFETIInstance *instance,
	esint size,
	double *values)
{
	std::vector<std::vector<double> > solution(1);
	solution[0] = std::vector<double>(values, values + size);
	(*instance)->data.solve(solution);
	memcpy(values, &solution[0][0], size);
	return 0;
}

int ESPRESODestroy(void *data)
{
	for (
		std::list<ESPRESOStructDoubleVector*>::iterator it = DataHolder::doubleVectors.begin();
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
		std::list<ESPRESOStructIntVector*>::iterator it = DataHolder::intVectors.begin();
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
		std::list<ESPRESOStructMap*>::iterator it = DataHolder::maps.begin();
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
		std::list<ESPRESOStructMatrix*>::iterator it = DataHolder::matrices.begin();
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
		std::list<ESPRESOStructFETIIntance*>::iterator it = DataHolder::instances.begin();
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

int ESPRESOFinalize()
{
	return 0;
}


