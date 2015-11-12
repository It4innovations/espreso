
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

struct ESPRESOStructIntVector {
	esint size;
	esint *values;
};

struct ESPRESOStructDoubleVector {
	esint size;
	double *values;
};

struct ESPRESOStructMap {
	esint size;
	esint *indices;
	double *values;
};


#ifdef __cplusplus
extern "C" {
#endif

typedef struct ESPRESOStructMat* ESPRESOMat;

typedef struct ESPRESOStructIntVector* ESPRESOIntVector;
typedef struct ESPRESOStructDoubleVector* ESPRESODoubleVector;
typedef struct ESPRESOStructMap* ESPRESOMap;

typedef struct ESPRESOStructFETIIntance* ESPRESOFETIInstance;


int ESPRESOInit(
	MPI_Comm communicator
);

int ESPRESOCreateMatrixElemental(
	esint n,
	esint nelt,
	esint *eltptr,
	esint *eltvar,
	double *values,
	ESPRESOMat *stiffnessMatrix
);

int ESPRESOCreateDoubleVector(
	esint size,
	double *values,
	ESPRESODoubleVector *vector
);

int ESPRESOCreateIntVector(
	esint size,
	esint *values,
	ESPRESOIntVector *vector
);

int ESPRESOCreateMap(
	esint size,
	esint *indices,
	double *values,
	ESPRESOMap *vector
);

int ESPRESOPrepareFETIInstance(
	ESPRESOMat *stiffnessMatrix,
	ESPRESODoubleVector *rhs,
	ESPRESOMap *dirichlet,
	ESPRESOIntVector *l2g,
	ESPRESOIntVector *neighbourRanks,
	ESPRESOFETIInstance *instance
);

int ESPRESOSolveFETI(
	ESPRESOFETIInstance *instance,
	esint size,
	double *values
);

int ESPRESODestroy(
	void *data
);

int ESPRESOFinalize();


#ifdef __cplusplus
}
#endif


#endif /* ESPRESO_API_H_ */
