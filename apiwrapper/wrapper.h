
#ifndef ESPRESO_WRAPPER_H_
#define ESPRESO_WRAPPER_H_

#include "../libespreso/espreso.h"
#include <iostream>

#include "esconfig.h"
#include "esassemblers.h"

struct ESPRESOMatData {
	SparseCSRMatrix<eslocal> &data;
};

struct DataHolder {
	static std::list<assembler::Assembler<assembler::API>*> assemblers;
	static MPI_Comm communicator;
};



#endif /* ESPRESO_WRAPPER_H_ */
