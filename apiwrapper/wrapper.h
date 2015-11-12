
#ifndef ESPRESO_WRAPPER_H_
#define ESPRESO_WRAPPER_H_

#include "../libespreso/espreso.h"
#include <iostream>

#include "esconfig.h"
#include "esassemblers.h"

struct ESPRESOStructMat {
	SparseCSRMatrix<eslocal> *data;
};

struct ESPRESOStructFETIIntance {
	assembler::API *data;
};

struct DataHolder {
	static std::list<ESPRESOStructMat> matrices;
	static std::list<ESPRESOFETIInstance> instances;
	static MPI_Comm communicator;
};



#endif /* ESPRESO_WRAPPER_H_ */
