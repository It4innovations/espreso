
#ifndef ESPRESO_WRAPPER_H_
#define ESPRESO_WRAPPER_H_

#include <iostream>

#include "esconfig.h"
#include "esassemblers.h"

struct FETI4IStructMatrix {
	FETI4IStructMatrix(): data(0, 0) { };

	SparseCSRMatrix<eslocal> data;
};

struct FETI4IStructIntance {
	FETI4IStructIntance(assembler::LinearElasticity<assembler::API2> data): data(data) { };

	assembler::LinearElasticity<assembler::API2> data;
};

struct DataHolder {
	static std::list<FETI4IStructMatrix*> matrices;
	static std::list<FETI4IStructIntance*> instances;
};



#endif /* ESPRESO_WRAPPER_H_ */
