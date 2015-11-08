
#ifndef ESPRESO_WRAPPER_H_
#define ESPRESO_WRAPPER_H_

#include "../libespreso/espreso.h"
#include <iostream>

#include "esassemblers.h"

struct DataHolder {
	static std::list<assembler::Assembler<assembler::API>*> assemblers;
	static MPI_Comm communicator;
	static int MPIrank;
	static int MPIsize;
};



#endif /* ESPRESO_WRAPPER_H_ */
