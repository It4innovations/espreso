
#include "wrapper.h"

using namespace assembler;

std::list<ESPRESOStructMat*> DataHolder::matrices;
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

	ESPRESOStructMat *holder = new ESPRESOStructMat(new SparseCSRMatrix<eslocal>(matrix));

	DataHolder::matrices.push_back(holder);
	*stiffnessMatrix = DataHolder::matrices.back();
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
//	std::list<ESPRESOStructMat*>::iterator it;
//	for (it = DataHolder::matrices.begin(); it != DataHolder::matrices.end(); ++it) {
//		if (*it == *stiffnessMatrix) {
//			std::cout << "MATRIX found\n";
//		}
//	}

	API api(*(*stiffnessMatrix)->data, **rhs, **dirichlet, **l2g, **neighbourRanks);
	ESPRESOStructFETIIntance *holder = new ESPRESOStructFETIIntance(new LinearElasticity<API>(api));

	holder->data->init();

	DataHolder::instances.push_back(holder);
	*instance = DataHolder::instances.back();
	return 0;
}

int ESPRESOSolveFETI(
	ESPRESOFETIInstance *instance,
	esint size,
	double *values)
{
//	std::list<ESPRESOStructFETIIntance*>::iterator it;
//	for (it = DataHolder::instances.begin(); it != DataHolder::instances.end(); ++it) {
//		if (*it == *instance) {
//			std::cout << "HURAA\n";
//		}
//	}

	std::vector<std::vector<double> > solution(1);
	solution[0] = std::vector<double>(values, values + size);
	(*instance)->data->solve(solution);
	memcpy(values, &solution[0][0], size);
	return 0;
}

int ESPRESODestroy(void *data)
{
	// TODO: implement me
	return 0;
}

int ESPRESOFinalize()
{
	std::list<ESPRESOStructMat*>::iterator it;
	for (it = DataHolder::matrices.begin(); it != DataHolder::matrices.end(); ++it) {
		delete *it;
	}
	DataHolder::matrices.clear();
	return 0;
}


