
#include "wrapper.h"

using namespace assembler;

std::list<ESPRESOStructMat> DataHolder::matrices;
std::list<ESPRESOFETIInstance> DataHolder::instances;
MPI_Comm DataHolder::communicator;

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
	ESPRESOMat *stiffnessMatrix)
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

	std::cout << matrix;
	return 0;
}

int ESPRESOCreateDoubleVector(
	esint size,
	double *values,
	ESPRESODoubleVector *vector)
{
	return 0;
}

int ESPRESOCreateIntVector(
	esint size,
	esint *values,
	ESPRESOIntVector *vector)
{
	return 0;
}

int ESPRESOCreateMap(
	esint size,
	esint *indices,
	double *values,
	ESPRESOMap *vector)
{
	return 0;
}

int ESPRESOPrepareFETIInstance(
	ESPRESOMat *stiffnessMatrix,
	ESPRESODoubleVector *rhs,
	ESPRESOMap *dirichlet,
	ESPRESOIntVector *l2g,
	ESPRESOIntVector *neighbourRanks,
	ESPRESOFETIInstance *instance)
{
	return 0;
}

int ESPRESOSolveFETI(
	ESPRESOFETIInstance *instance,
	esint size,
	double *values)
{
	return 0;
}

int ESPRESODestroy(void *data)
{
	// TODO: implement me
	return 0;
}

int ESPRESOFinalize()
{
//	std::list<Assembler<API>*>::iterator it;
//	for (it = DataHolder::assemblers.begin(); it != DataHolder::assemblers.end(); ++it) {
//		delete *it;
//	}
//	DataHolder::assemblers.clear();
	return 0;
}


