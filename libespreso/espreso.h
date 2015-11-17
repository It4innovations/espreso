
#ifndef ESPRESO_API_H_
#define ESPRESO_API_H_

#include "mpi.h"

/*-----------------------------------------------------------------------------
 Set data-types used in ESPRESO

 Possible values for FETI4I_INDICES_WIDTH:
   32: use 32 bit signed integer
   64: use 64 bit signed integer

 Possible values for FETI4I_REAL_WIDTH:
   64: ESPRESO supports only 64 bit real values
------------------------------------------------------------------------------*/
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


#ifdef __cplusplus
extern "C" {
#endif

/*-----------------------------------------------------------------------------
 Definitions of internal structures used in ESPRESO
------------------------------------------------------------------------------*/
typedef struct FETI4IStructMatrix* FETI4IMatrix;
typedef struct FETI4IStructIntVector* FETI4IIntVector;
typedef struct FETI4IStructDoubleVector* FETI4IDoubleVector;
typedef struct FETI4IStructMap* FETI4IMap;

typedef struct FETI4IStructFETIIntance* FETI4IFETIInstance;


/*-----------------------------------------------------------------------------
 Functions for manipulating with ESPRESO internal structures
------------------------------------------------------------------------------*/
int FETI4ICreateMatrixElemental(
	FETI4IInt n,
	FETI4IInt nelt,
	FETI4IInt *eltptr,
	FETI4IInt *eltvar,
	FETI4IReal *values,
	FETI4IMatrix *stiffnessMatrix
);

int FETI4ICreateDoubleVector(
	FETI4IInt size,
	FETI4IReal *values,
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
	FETI4IReal *values,
	FETI4IMap *vector
);

int FETI4IUpdateDoubleVector(
	FETI4IInt size,
	FETI4IReal *values,
	FETI4IDoubleVector *vector
);

int FETI4IUpdateIntVector(
	FETI4IInt size,
	FETI4IInt *values,
	FETI4IIntVector *vector
);

int FETI4IUpdateMap(
	FETI4IInt size,
	FETI4IInt *indices,
	FETI4IReal *values,
	FETI4IMap *vector
);


/*-----------------------------------------------------------------------------
 Functions for creating an instance and solve it
------------------------------------------------------------------------------*/

int FETI4ISolveFETI(
	FETI4IInt *settings,
	FETI4IFETIInstance *instance,
	FETI4IInt size,
	FETI4IReal *values
);

int FETI4IPrepareFETIInstance(
	FETI4IInt *settings,
	FETI4IMatrix *stiffnessMatrix,
	FETI4IDoubleVector *rhs,
	FETI4IMap *dirichlet,
	FETI4IIntVector *l2g,
	FETI4IIntVector *neighbourRanks,
	FETI4IFETIInstance *instance
);

int FETI4IUpdateStiffnessMatrix(
	FETI4IMatrix *stiffnessMatrix,
	FETI4IFETIInstance *instance
);

int FETI4IUpdateRhs(
	FETI4IDoubleVector *rhs,
	FETI4IFETIInstance *instance
);

int FETI4IUpdateDirichlet(
	FETI4IMap *dirichlet,
	FETI4IIntVector *l2g,
	FETI4IIntVector *neighbourRanks,
	FETI4IFETIInstance *instance
);


/*-----------------------------------------------------------------------------
 Destroy an arbitrary internal structure
------------------------------------------------------------------------------*/

int FETI4IDestroy(
	void *data
);



////// Avoid Init and Finalize to make more user-friendly approach?
int FETI4IInit(
	MPI_Comm communicator
);
int FETI4IFinalize();
///////////////////////////////

#ifdef __cplusplus
}
#endif


#endif /* ESPRESO_API_H_ */
