
#ifndef ESPRESO_WRAPPER_H_
#define ESPRESO_WRAPPER_H_

#include "../libespreso/espreso.h"
#include <iostream>

#include "esconfig.h"
#include "esassemblers.h"

struct ESPRESOStructMat {
	ESPRESOStructMat(SparseCSRMatrix<eslocal> *data): data(data) {};
	~ESPRESOStructMat() { delete data; }

	SparseCSRMatrix<eslocal> *data;
};

struct ESPRESOStructFETIIntance {
	ESPRESOStructFETIIntance(assembler::LinearElasticity<assembler::API> *data): data(data) {};
	~ESPRESOStructFETIIntance() { delete data; }

	assembler::LinearElasticity<assembler::API> *data;
};

struct DataHolder {
	static std::list<ESPRESOStructMat*> matrices;
	static std::list<ESPRESOStructFETIIntance*> instances;
	static MPI_Comm communicator;
};



#endif /* ESPRESO_WRAPPER_H_ */
