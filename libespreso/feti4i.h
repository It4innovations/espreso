
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
/// Data holder for a stiffness matrix
typedef struct FETI4IStructMatrix* FETI4IMatrix;
/// Data holder for an instance of a problem solvable by the ESPRESO solver
typedef struct FETI4IStructInstance* FETI4IInstance;

/*-----------------------------------------------------------------------------
 Functions for manipulating with FETI4I internal structures
------------------------------------------------------------------------------*/

/// Create the data holder for a stiffness matrix
/**
 * Create an ESPRESO internal representation of a stiffness matrix.
 * Created matrix is an empty matrix prepared for later adding of element matrices.
 *
 * @param matrix an address of a data holder
 * @param indexBase a value of the first matrix index
 */
void FETI4ICreateStiffnessMatrix(
		FETI4IMatrix 	*matrix,
		FETI4IInt		indexBase
);

/// Adds an element to a stiffness matrix
/**
 * Adds a new element to a stiffness matrix.
 * The method assumes symmetric stiffness element matrix.
 *
 * @param matrix a data holder of a stiffness matrix
 * @param size the number of rows of square matrix
 * @param indices an array of local DOFs of an element matrix [size]
 * @param values an array of values in row-major-order [size x size]
 */
void FETI4IAddElement(
		FETI4IMatrix 	matrix,
		FETI4IInt 		size,
		FETI4IInt* 		indices,
		FETI4IReal* 	values
);

/*-----------------------------------------------------------------------------
 Functions for creating an instance and solve it
------------------------------------------------------------------------------*/

/// Create the data holder for a problem solvable in the ESPRESO solver
/**
 * Create an ESPRESO internal representation of a problem.
 * Created instance can be directly passed to the ESPRESO solver
 *
 * @param instance an address of a data holder of an ESPRESO instance
 * @param matrix a data holder of a stiffness matrix
 * @param size size of a right-hand side and a mapping of local DOFs to global DOFs
 * @param rhs right-hand side
 * @param l2g a mapping of local DOFs to global DOFs
 * @param neighbours_size a number of neighbours MPI ranks
 * @param neighbours an array of neighbours MPI ranks
 * @param dirichlet_size a number of DOFs with dirichlet condition
 * @param dirichlet_indices an array of DOFs (local)
 * @param dirichlet_values an array of dirichlet value
 */
void FETI4ICreateInstance(
		FETI4IInstance 	*instance,
		FETI4IMatrix 	matrix,
		FETI4IInt 		size,
		FETI4IReal* 	rhs,
		FETI4IInt* 		l2g,
		FETI4IMPIInt 	neighbours_size,
		FETI4IMPIInt*	neighbours,
		FETI4IInt 		dirichlet_size,
		FETI4IInt* 		dirichlet_indices,
		FETI4IReal* 	dirichlet_values
);


/// Solve an instance by the ESPRESO solver
/**
 * The method gets an instance and returns the solution to a pre-allocated array
 *
 * @param instance a data holder of an ESPRESO instance
 * @param solution_size a solution size
 * @param solution an array where is save the solution
 */
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

/// Destroy a data holder
/**
 * Destroy all data associated to a data holder
 *
 * @param ptr an address to a data holder
 */
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

