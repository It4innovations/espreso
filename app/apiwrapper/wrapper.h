
#ifndef ESPRESO_WRAPPER_H_
#define ESPRESO_WRAPPER_H_

#include "../libespreso/espreso.h"
#include <iostream>

#include "esconfig.h"
#include "esassemblers.h"

struct ESPRESOStructMatrix {
	ESPRESOStructMatrix(): data(0, 0) { };

	SparseCSRMatrix<eslocal> data;
};

struct ESPRESOStructFETIIntance {
	ESPRESOStructFETIIntance(assembler::LinearElasticity<assembler::API> data): data(data) { };

	assembler::LinearElasticity<assembler::API> data;
};

struct DataHolder {
	static std::list<ESPRESOStructDoubleVector*> doubleVectors;
	static std::list<ESPRESOStructIntVector*> intVectors;
	static std::list<ESPRESOStructMap*> maps;
	static std::list<ESPRESOStructMatrix*> matrices;
	static std::list<ESPRESOStructFETIIntance*> instances;
	static MPI_Comm communicator;
};



#endif /* ESPRESO_WRAPPER_H_ */
