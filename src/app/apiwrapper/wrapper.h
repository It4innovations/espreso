
#ifndef ESPRESO_WRAPPER_H_
#define ESPRESO_WRAPPER_H_

#include <iostream>
#include <list>

#include "esconfig.h"
#include "esassembler.h"
#include "esinput.h"
#include "esbasis.h"

struct FETI4IStructMatrix {
	FETI4IStructMatrix(eslocal offset): offset(offset) {};

	std::vector<std::vector<eslocal> > eIndices;
	std::vector<std::vector<double> > eMatrices;

	eslocal offset;
};

struct FETI4IStructInstance {
	FETI4IStructInstance(FETI4IStructMatrix &matrix)
	: instance(NULL), mesh(matrix.eMatrices) {};
	~FETI4IStructInstance() { if (instance != NULL) { delete instance; } }

	espreso::Instance *instance;
	espreso::APIMesh mesh;
};

namespace espreso {

struct DataHolder {
	static std::list<FETI4IStructMatrix*> matrices;
	static std::list<FETI4IStructInstance*> instances;
	static TimeEval timeStatistics;
};

}


#endif /* ESPRESO_WRAPPER_H_ */
