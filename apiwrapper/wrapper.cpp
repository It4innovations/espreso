
#include "wrapper.h"

using namespace assembler;

std::list<Assembler<API>*> DataHolder::assemblers;
MPI_Comm DataHolder::communicator;

int ESPRESOInit(MPI_Comm communicator)
{
	DataHolder::communicator = communicator;
	MPI_Comm_rank(communicator, &esconfig::MPIrank);
	MPI_Comm_size(communicator, &esconfig::MPIsize);
	return 0;
}

int ESPRESOCreateStiffnessMatrix(
	esint n,
	esint nelms,
	esint *eltptr,
	esint *eltvar,
	double *values,
	ESPRESOMat *stiffnessMatrix)
{
	return 0;
}

int ESPRESOSolveFETI(
	ESPRESOMat *stiffnessMatrix,
	ESPRESODoubleVector *rhs,
	ESPRESOMap *dirichlet,
	ESPRESOIntVector *l2g,
	ESPRESOIntVector *neighbourRanks,
	ESPRESODoubleVector *solution)
{
	return 0;
}

int ESPRESOFree(void *data)
{
	// TODO: implement me
	return 0;
}

int ESPRESOFinalize()
{
	std::list<Assembler<API>*>::iterator it;
	for (it = DataHolder::assemblers.begin(); it != DataHolder::assemblers.end(); ++it) {
		delete *it;
	}
	DataHolder::assemblers.clear();
	return 0;
}


