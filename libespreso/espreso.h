
#ifndef ESPRESO_API_H_
#define ESPRESO_API_H_

#include "mpi.h"

#define FETI4I_INDICES_WIDTH 32
#define FETI4I_REAL_WIDTH 64


#if FETI4I_INDICES_WIDTH == 32
	typedef int FETI4IInt;
#elif FETI4I_INDICES_WIDTH == 64
	typedef long esint;
#else
	#error "Incorrect user-supplied value of FETI4I_INDICES_WIDTH"
#endif

#if FETI4I_REAL_WIDTH == 64
	typedef long FETI4IReal;
#else
	#error "Incorrect user-supplied value of FETI4I_REAL_WIDTH"
#endif


struct FETI4IStructIntVector {
	FETI4IInt size;
	FETI4IInt *values;
};

struct FETI4IStructDoubleVector {
	FETI4IInt size;
	double *values;
};

struct FETI4IStructMap {
	FETI4IInt size;
	FETI4IInt *indices;
	double *values;
};

#ifdef __cplusplus
extern "C" {
#endif

typedef struct FETI4IStructMatrix* FETI4IMatrix;
typedef struct FETI4IStructIntVector* FETI4IIntVector;
typedef struct FETI4IStructDoubleVector* FETI4IDoubleVector;
typedef struct FETI4IStructMap* FETI4IMap;

typedef struct FETI4IStructFETIIntance* FETI4IFETIInstance;

////// Avoid init and finalize?
int FETI4IInit(
	MPI_Comm communicator
);
int FETI4IFinalize();
///////////////////////////////

int FETI4ICreateMatrixElemental(
	FETI4IInt n,
	FETI4IInt nelt,
	FETI4IInt *eltptr,
	FETI4IInt *eltvar,
	double *values,
	FETI4IMatrix *stiffnessMatrix
);

int FETI4ICreateDoubleVector(
	FETI4IInt size,
	double *values,
	FETI4IDoubleVector *vector
);

int FETI4ICreateIntVector(
	FETI4IInt size,
	FETI4IInt *values,
	FETI4IIntVector *vector
);

int FETI4ICreateMap(
	FETI4IInt size,
	FETI4IInt *indices,
	double *values,
	FETI4IMap *vector
);

int FETI4IPrepareFETIInstance(
	FETI4IMatrix *stiffnessMatrix,
	FETI4IDoubleVector *rhs,
	FETI4IMap *dirichlet,
	FETI4IIntVector *l2g,
	FETI4IIntVector *neighbourRanks,
	FETI4IFETIInstance *instance
);

int FETI4ISolveFETI(
	FETI4IFETIInstance *instance,
	FETI4IInt size,
	double *values
);

int FETI4IDestroy(
	void *data
);


#ifdef __cplusplus
}
#endif


#endif /* ESPRESO_API_H_ */
