
#ifndef FETI4I_H_
#define FETI4I_H_

#include "mpi.h"

/*-----------------------------------------------------------------------------
 Set data-types used in FETI4I

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
	typedef double FETI4IReal;
#else
	#error "Incorrect user-supplied value of FETI4I_REAL_WIDTH"
#endif


#ifdef __cplusplus
extern "C" {
#endif

/*-----------------------------------------------------------------------------
 Definitions of internal structures used in FETI4I
------------------------------------------------------------------------------*/
typedef struct FETI4IStructMatrix* FETI4IMatrix;
typedef struct FETI4IStructIntance* FETI4IInstance;


/*-----------------------------------------------------------------------------
 Functions for manipulating with FETI4I internal structures
------------------------------------------------------------------------------*/
int FETI4ICreateStiffnessMatrix(
		FETI4IMatrix 	*stiffnessMatrix,
		FETI4IInt 		n,
		FETI4IInt 		nelt,
		FETI4IInt* 		eltptr,
		FETI4IInt* 		eltvar,
		FETI4IReal* 	values
);

/*-----------------------------------------------------------------------------
 Functions for creating an instance and solve it
------------------------------------------------------------------------------*/

int FETI4ISolve(
		FETI4IInstance 	instance,
		FETI4IInt 		solution_size,
		FETI4IReal*		solution
);

int FETI4ICreateInstance(
		FETI4IInstance 	*instance,
		FETI4IInt* 		settings,	// Currently only NULL is supported
		FETI4IMatrix 	stiffnessMatrix,
		FETI4IInt 		rhs_size,
		FETI4IReal* 	rhs,
		FETI4IInt 		dirichlet_size,
		FETI4IInt* 		dirichlet_indices,
		FETI4IReal* 	dirichlet_values,
		FETI4IInt 		l2g_size,
		FETI4IInt* 		l2g,
		FETI4IInt 		neighbours_size,
		FETI4IInt* 		neighbours,
		MPI_Comm 		communicator	// Currently only MPI_COMM_WORLD is supported
);

/*-----------------------------------------------------------------------------
 Functions for updating a created instance
------------------------------------------------------------------------------*/

int FETI4IUpdateStiffnessMatrix(
		FETI4IInstance 	instance,
		FETI4IMatrix 	stiffnessMatrix
);

int FETI4IUpdateRhs(
		FETI4IInstance 	instance,
		FETI4IInt 		rhs_size,
		FETI4IReal* 	rhs_values
);

int FETI4IUpdateDirichlet(
		FETI4IInstance 	instance,
		FETI4IInt 		dirichlet_size,
		FETI4IInt* 		dirichlet_indices,
		FETI4IReal* 	dirichlet_values,
		FETI4IInt 		l2g_size,
		FETI4IReal* 	l2g_values,
		FETI4IInt 		neighbour_size,
		FETI4IReal* 	neighbour_values
);


/*-----------------------------------------------------------------------------
 Destroy an arbitrary internal structure
------------------------------------------------------------------------------*/

int FETI4IDestroy(
		void* 			data
);

#ifdef __cplusplus
}
#endif


#endif /* FETI4I_H_ */
