
#ifndef ESPRESO_API_H_
#define ESPRESO_API_H_

#include "mpi.h"

#define ESPRESO_INDICES_WIDTH 32


#if ESPRESO_INDICES_WIDTH == 32
	typedef int esint;
#elif ESPRESO_INDICES_WIDTH == 64
	typedef long esint;
#else
	#error "Incorrect user-supplied value of ESPRESO_INDICES_WIDTH"
#endif


struct LocalStiffnessMatrices {
	esint size;
	esint *array;
};

struct DoubleVector {
	esint size;
	double *values;
};

struct IntVector {
	esint size;
	esint *values;
};

#ifdef __cplusplus
extern "C" {
#endif

typedef struct LocalStiffnessMatrices LocalStiffnessMatrices;
typedef struct DoubleVector DoubleVector;
typedef void* ESPRESOHandler;

int ESPRESO_Init(
	MPI_Comm communicator
);

int FETI_PrepareElasticity(
	LocalStiffnessMatrices *Ke,
	DoubleVector *rhs,
	IntVector *l2g,
	ESPRESOHandler *handler
);

int Solve(
	ESPRESOHandler *handler,
	DoubleVector *solution
);

int ESPRESO_Finalize();


#ifdef __cplusplus
}
#endif


#endif /* ESPRESO_API_H_ */
