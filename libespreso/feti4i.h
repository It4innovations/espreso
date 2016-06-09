
#ifndef FETI4I_H_
#define FETI4I_H_

#include "mpi.h"

/*-----------------------------------------------------------------------------
 Set data-types used in FETI4I

 Possible values for FETI4I_INT_WIDTH:
   32: use 32 bit signed integer
   64: use 64 bit signed integer

 Possible values for FETI4I_REAL_WIDTH:
   64: ESPRESO supports only 64 bit real values
------------------------------------------------------------------------------*/
#ifndef FETI4I_INT_WIDTH
#define FETI4I_INT_WIDTH 32
#endif

#ifndef FETI4I_REAL_WIDTH
#define FETI4I_REAL_WIDTH 64
#endif

#if FETI4I_INT_WIDTH == 32
	typedef int FETI4IInt;
#elif FETI4I_INT_WIDTH == 64
	typedef long FETI4IInt;
#else
#error "Incorrect user-supplied value of FETI4I_INT_WIDTH"
#endif

/* MPI integer (e.g. rank) are always 32-bit */
	typedef int FETI4IMPIInt;

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
typedef struct FETI4IStructInstance* FETI4IInstance;

/*-----------------------------------------------------------------------------
 Functions for manipulating with FETI4I internal structures
------------------------------------------------------------------------------*/

void FETI4ICreateStiffnessMatrix(
		FETI4IMatrix 	*matrix,           // pointer to a created matrix
		FETI4IInt		indexBase          // value of the first index
);

void FETI4IAddElement(
		FETI4IMatrix 	matrix,            // pointer to a matrix
		FETI4IInt 		size,              // number of rows of square element matrix
		FETI4IInt* 		indices,           // vector of indices (size)
		FETI4IReal* 	values             // vector of values (size x size)
);

/*-----------------------------------------------------------------------------
 Functions for creating an instance and solve it
------------------------------------------------------------------------------*/

void FETI4ICreateInstance(
		FETI4IInstance 	*instance,         // pointer to a created instance
		FETI4IMatrix 	matrix,            // stiffness matrix
		FETI4IInt 		size,              // size of rhs and l2g
		FETI4IReal* 	rhs,               // right-hand side
		FETI4IInt* 		l2g,               // local to global mapping
		FETI4IMPIInt 	neighbours_size,   // number of neighbours
		FETI4IMPIInt*	neighbours,        // vector of neighbours
		FETI4IInt 		dirichlet_size,    // number of DOFs with dirichlet
		FETI4IInt* 		dirichlet_indices, // vector of dirichlet indices in local indexing
		FETI4IReal* 	dirichlet_values   // vector of dirichlet values
);

void FETI4ISolve(
		FETI4IInstance 	instance,          // pointer to an instance
		FETI4IInt 		solution_size,     // size of solution vector
		FETI4IReal*		solution           // solution
);


/*-----------------------------------------------------------------------------
 Functions for updating a created instance
------------------------------------------------------------------------------*/

void FETI4IUpdateStiffnessMatrix(
		FETI4IInstance 	instance,
		FETI4IMatrix 	stiffnessMatrix
);

void FETI4IUpdateRhs(
		FETI4IInstance 	instance,
		FETI4IInt 		rhs_size,
		FETI4IReal* 	rhs_values
);

//TODO VH: neighbours perhaps should not be passed here and they are FETI4IMPIInt
void FETI4IUpdateDirichlet(
		FETI4IInstance 	instance,
		FETI4IInt 		dirichlet_size,
		FETI4IInt* 		dirichlet_indices,
		FETI4IReal* 	dirichlet_values,
		FETI4IReal* 	l2g,
		FETI4IInt 		neighbour_size,
		FETI4IReal* 	neighbour
);


/*-----------------------------------------------------------------------------
 Destroy an arbitrary internal structure
------------------------------------------------------------------------------*/

void FETI4IDestroy(
		void* 			ptr
);

/*-----------------------------------------------------------------------------
 API testing functions
------------------------------------------------------------------------------*/

void FETI4ITest();

void TEST4IGetElementsInfo(
		FETI4IInt		*elements,
		FETI4IInt		*elementSize
);

void TEST4IGetElement(
		FETI4IInt		index,
		FETI4IInt*		*indices,
		FETI4IReal*		*values
);

void TEST4IGetInstanceInfo(
		FETI4IInt		*rhs_size,
		FETI4IInt		*dirichlet_size,
		FETI4IInt		*neighbours_size
);

void TEST4IGetInstance(
		FETI4IReal*		*rhs,
		FETI4IInt*		*l2g,
		FETI4IInt*		*dirichlet_indices,
		FETI4IReal*		*dirichlet_values,
		FETI4IInt*		*neighbours
);


#ifdef __cplusplus
}
#endif


#endif /* FETI4I_H_ */

