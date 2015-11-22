
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
#ifndef FETI4I_INDICES_WIDTH
#define FETI4I_INDICES_WIDTH 32
#endif

#define FETI4I_REAL_WIDTH 64

#if FETI4I_INDICES_WIDTH == 32
	typedef int FETI4IInt;
#elif FETI4I_INDICES_WIDTH == 64
	typedef long FETI4IInt;
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
typedef struct FETI4IStructRealVector* FETI4IRealVector;

typedef struct FETI4IStructIntance* FETI4IInstance;


/*-----------------------------------------------------------------------------
 Functions for manipulating with ESPRESO internal structures
------------------------------------------------------------------------------*/
int FETI4ICreateMatrixElemental(
		FETI4IMatrix *stiffnessMatrix,
		FETI4IInt n,
		FETI4IInt nelt,
		FETI4IInt* eltptr,
		FETI4IInt* eltvar,
		FETI4IReal* values
);

int FETI4ICreateDoubleVector(
		FETI4IRealVector *vector,
		FETI4IInt size,
		FETI4IReal* values
);

int FETI4ICreateIntVector(
		FETI4IIntVector *vector,
		FETI4IInt size,
		FETI4IInt *values
);

int FETI4IUpdateRealVector(
		FETI4IRealVector vector,
		FETI4IInt size,
		FETI4IReal* values
);

int FETI4IUpdateIntVector(
		FETI4IIntVector vector,
		FETI4IInt size,
		FETI4IInt* values
);


/*-----------------------------------------------------------------------------
 Functions for creating an instance and solve it
------------------------------------------------------------------------------*/

int FETI4ISolve(
	FETI4IInstance instance,
	FETI4IInt solution_size,
	FETI4IReal* solution
);

int FETI4ICreateInstance(
		FETI4IInstance *instance,
		FETI4IInt* settings,
		FETI4IMatrix stiffnessMatrix,
		FETI4IRealVector rhs,
		FETI4IIntVector dirichlet_indices,
		FETI4IRealVector dirichlet_values,
		FETI4IIntVector l2g,
		FETI4IIntVector neighbourRanks,
		MPI_Comm communicator
);

int FETI4IUpdateStiffnessMatrix(
		FETI4IInstance instance,
		FETI4IMatrix stiffnessMatrix
);

int FETI4IUpdateRhs(
		FETI4IInstance instance,
		FETI4IRealVector rhs
);

int FETI4IUpdateDirichlet(
		FETI4IInstance instance,
		FETI4IIntVector dirichlet_instances,
		FETI4IRealVector dirichlet_values,
		FETI4IIntVector l2g,
		FETI4IIntVector neighbourRanks
);


/*-----------------------------------------------------------------------------
 Destroy an arbitrary internal structure
------------------------------------------------------------------------------*/

int FETI4IDestroy(
	void *data
);

#ifdef __cplusplus
}
#endif


#endif /* ESPRESO_API_H_ */
