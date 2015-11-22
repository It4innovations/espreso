
#ifndef ESPRESO_WRAPPER_H_
#define ESPRESO_WRAPPER_H_

#include <iostream>

#include "esconfig.h"
#include "esassemblers.h"

struct FETI4IStructRealVector {
	std::vector<double> data;
};

struct FETI4IStructIntVector {
	std::vector<eslocal> data;
};

struct FETI4IStructMatrix {
	FETI4IStructMatrix(): data(0, 0) { };

	SparseCSRMatrix<eslocal> data;
};

struct FETI4IStructIntance {
	FETI4IStructIntance(assembler::LinearElasticity<assembler::API> data): data(data) { };

	assembler::LinearElasticity<assembler::API> data;
};

struct DataHolder {
	static std::list<FETI4IStructRealVector*> doubleVectors;
	static std::list<FETI4IStructIntVector*> intVectors;
	static std::list<FETI4IStructMatrix*> matrices;
	static std::list<FETI4IStructIntance*> instances;
	static MPI_Comm communicator;
};



#endif /* ESPRESO_WRAPPER_H_ */
