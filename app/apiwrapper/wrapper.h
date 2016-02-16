
#ifndef ESPRESO_WRAPPER_H_
#define ESPRESO_WRAPPER_H_

#include <iostream>
#include <list>

#include "esconfig.h"
#include "esassemblers.h"
#include "esinput.h"

struct FETI4IStructMatrix {
	FETI4IStructMatrix(eslocal offset): offset(offset) { };

	std::vector<SparseVVPMatrix<eslocal> > K;

	std::vector<std::vector<eslocal> > eIndices;
	std::vector<std::vector<double> > eMatrix;

	eslocal offset;
};

struct FETI4IStructInstance {
	FETI4IStructInstance(assembler::LinearElasticity<assembler::API> data): data(data) { };

	assembler::LinearElasticity<assembler::API> data;
	std::vector<SparseCSRMatrix<eslocal> > K;
};

struct DataHolder {
	static std::list<FETI4IStructMatrix*> matrices;
	static std::list<FETI4IStructInstance*> instances;
};



#endif /* ESPRESO_WRAPPER_H_ */
