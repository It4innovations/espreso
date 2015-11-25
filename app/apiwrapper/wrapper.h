
#ifndef ESPRESO_WRAPPER_H_
#define ESPRESO_WRAPPER_H_

#include <iostream>

#include "esconfig.h"
#include "esassemblers.h"
#include "esinput.h"

struct FETI4IStructMatrix {
	FETI4IStructMatrix(eslocal offset): data(0, 0), offset(offset) { };

	SparseVVPMatrix<eslocal> data;
	eslocal offset;
};

struct FETI4IStructInstance {
	FETI4IStructInstance(assembler::LinearElasticity<assembler::API2> data): data(data), K(0, 0) { };

	assembler::LinearElasticity<assembler::API2> data;
	SparseCSRMatrix<eslocal> K;
};

struct DataHolder {
	static std::list<FETI4IStructMatrix*> matrices;
	static std::list<FETI4IStructInstance*> instances;
};



#endif /* ESPRESO_WRAPPER_H_ */
