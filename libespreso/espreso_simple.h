
#ifndef ESPRESO_API_SIMPLE_H_
#define ESPRESO_API_SIMPLE_H_

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

/*-----------------------------------------------------------------------------
 Functions for creating an instance and solve it
------------------------------------------------------------------------------*/

int FETI4ISolveFETI(
	FETI4IInt *settings,
	FETI4IFETIInstance *instance,
	FETI4IInt solution_size,
	FETI4IReal *solution_values
);

int FETI4IPrepareFETIInstance(
	FETI4IInt *settings,
	FETI4IMatrix *stiffnessMatrix,
	FETI4IInt rhs_size,
	FETI4IReal *rhs_values,
	FETI4IInt dirichlet_size,
	FETI4IInt *dirichlet_indices,
	FETI4IReal *dirichlet_values,
	FETI4IInt l2g_size,
	FETI4IReal *l2g_values,
	FETI4IInt neighbour_size,
	FETI4IReal *neighbour_values,
	MPI_Comm communicator,
	FETI4IFETIInstance *instance
);

int FETI4IUpdateStiffnessMatrix(
	FETI4IMatrix *stiffnessMatrix,
	FETI4IFETIInstance *instance
);

int FETI4IUpdateRhs(
	FETI4IInt rhs_size,
	FETI4IReal *rhs_values,
	FETI4IFETIInstance *instance
);

int FETI4IUpdateDirichlet(
	FETI4IInt dirichlet_size,
	FETI4IInt *dirichlet_indices,
	FETI4IReal *dirichlet_values,
	FETI4IInt l2g_size,
	FETI4IReal *l2g_values,
	FETI4IInt neighbour_size,
	FETI4IReal *neighbour_values,
	FETI4IFETIInstance *instance
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


#endif /* ESPRESO_API_SIMPLE_H_ */
