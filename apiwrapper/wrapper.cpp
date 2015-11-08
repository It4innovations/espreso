
#include "wrapper.h"

using namespace assembler;

std::list<Assembler<API>*> DataHolder::assemblers;
MPI_Comm DataHolder::communicator;
int DataHolder::MPIrank;
int DataHolder::MPIsize;


int ESPRESO_Init(MPI_Comm communicator)
{
	DataHolder::communicator = communicator;
	MPI_Comm_rank(communicator, &DataHolder::MPIrank);
	MPI_Comm_size(communicator, &DataHolder::MPIsize);
	return 0;
}

int ESPRESO_Finalize()
{
	std::list<Assembler<API>*>::iterator it;
	for (it = DataHolder::assemblers.begin(); it != DataHolder::assemblers.end(); ++it) {
		delete *it;
	}
	DataHolder::assemblers.clear();
	return 0;
}

int FETI_PrepareElasticity(
	LocalStiffnessMatrices *Ke,
	DoubleVector *rhs,
	ESPRESOHandler *handler)
{
	API api(*Ke, *rhs);
	DataHolder::assemblers.push_back(new LinearElasticity<API>(api));
	*handler = DataHolder::assemblers.back();
	return 0;
}

int Solve(
	ESPRESOHandler *handler,
	DoubleVector *solution)
{
	std::list<Assembler<API>*>::iterator it;
	for (it = DataHolder::assemblers.begin(); it != DataHolder::assemblers.end(); ++it) {
		if (*it == *handler) {
			// it->solve(...)
		}
	}
	return 0;
}



