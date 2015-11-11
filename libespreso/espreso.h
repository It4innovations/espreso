
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

struct ESPRESOIntVector {
	esint size;
	esint *values;
};

struct ESPRESODoubleVector {
	esint size;
	double *values;
};

struct ESPRESOMap {
	esint size;
	esint *indices;
	double *values;
};


#ifdef __cplusplus
extern "C" {
#endif

typedef struct ESPRESOMatData* ESPRESOMat;

typedef struct ESPRESOIntVector ESPRESOIntVector;
typedef struct ESPRESODoubleVector ESPRESODoubleVector;
typedef struct ESPRESOMap ESPRESOMap;


int ESPRESOInit(
	MPI_Comm communicator
);

int ESPRESOCreateStiffnessMatrix(
	esint n,
	esint nelms,
	esint *eltptr,
	esint *eltvar,
	double *values,
	ESPRESOMat *stiffnessMatrix
);

int ESPRESOSolveFETI(
	ESPRESOMat *stiffnessMatrix,
	ESPRESODoubleVector *rhs,
	ESPRESOMap *dirichlet,
	ESPRESOIntVector *l2g,
	ESPRESOIntVector *neighbourRanks,
	ESPRESODoubleVector *solution
);

int ESPRESOFree(
	void *data
);

int ESPRESOFinalize();


#ifdef __cplusplus
}
#endif


#endif /* ESPRESO_API_H_ */
