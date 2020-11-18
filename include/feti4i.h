
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

/// Initialize the FETI4I library
/**
 * ESPRESO needs to be initialized before calling any other function.
 * At the end of the program FETI4IFinalize should be called.
 *
 * @param comm an MPI communication
 * @param verbosity a verbosity level
 */
void FETI4IInit(
		MPI_Comm 		comm,
		FETI4IInt		verbosity
);

/// Finalize the FETI4I library
/**
 * This function destroy all internal parameters.
 */
void FETI4IFinalize();

/// Fill integer options to the default values
/**
 * ESPRESO is controlled by options passed to FETI4ICreateInstance method.
 * This method set all integer options to its default values.
 *
 * @param options an array of size FETI4I_INTEGER_OPTIONS_SIZE
 */
void FETI4ISetDefaultIntegerOptions(
		FETI4IInt*		options
);

/// Fill integer options to the default values
/**
 * ESPRESO is controlled by options passed to FETI4ICreateInstance method.
 * This method set all real options to its default values.
 *
 * @param options an array of size FETI4I_REAL_OPTIONS_SIZE
 */
void FETI4ISetDefaultRealOptions(
		FETI4IReal*		options
);

/// Create the data holder for a stiffness matrix
/**
 * Create an ESPRESO internal representation of a stiffness matrix.
 * Created matrix is an empty matrix prepared for later adding of element matrices.
 *
 * @param matrix an address of a data holder
 * @param nodesSize a number of local nodes
 * @param l2g a mapping of local nodes (offsets) to global nodes (ids)
 * @param indexBase a value of the first matrix index
 * @param type a type of the matrix
 * @param dofsPerNode a number of DOFs per each nodes
 * @param dofsOrdering ordering of DOFs within each element (grouped by DOFs [xx.., yy.., zz..] or grouped by nodes [xyz, xyz, ...])
 */
void FETI4ICreateStiffnessMatrix(
		FETI4IMatrix	*matrix,
		FETI4IInt		nodesSize,
		FETI4IInt*		l2g,
		FETI4IInt		indexBase,
		FETI4IInt		type,
		FETI4IInt		dofsPerNode,
		FETI4IInt		dofsOrdering
);

/// Adds an element to a stiffness matrix
/**
 * Adds a new element to a stiffness matrix.
 * It is possible to call this method in parallel.
 *
 * @param matrix a data holder of a stiffness matrix
 * @param type a type of the element (currently only FETI4I_ETYPE_VOLUME is supported)
 * @param size the number of nodes of the added element
 * @param nodes an array of local nodes of the element matrix
 * @param stiffness stiffness values in row-major-order [(DOFs * size) x (DOFs * size)]
 */
void FETI4IAddElement(
		FETI4IMatrix	matrix,
		FETI4IInt		type,
		FETI4IInt		size,
		FETI4IInt*		nodes,
		FETI4IReal*		stiffness
);

/*-----------------------------------------------------------------------------
 Functions for creating an instance and solve it
------------------------------------------------------------------------------*/

/// Create the data holder for a problem solvable in the ESPRESO solver
/**
 * Create an ESPRESO internal representation of a problem.
 * Created instance can be directly passed to the ESPRESO solver or edited
 * by appropriate API methods.
 *
 * @param instance an address of a data holder of an ESPRESO instance
 * @param matrix a data holder of a stiffness matrix
 * @param neighbours_size a number of neighbours MPI ranks
 * @param neighbours an array of neighbours MPI ranks
 * @param dirichlet_size a number of DOFs with dirichlet condition
 * @param dirichlet_indices an array of DOFs (local)
 * @param dirichlet_values an array of dirichlet value
 * @param integer_options an array of integer parameters - see enum 'FETI4IIntegerOptions'
 * @param real_options an array of real parameters - see enum 'FETI4IRealOptions'
 */
void FETI4ICreateInstance(
		FETI4IInstance	*instance,
		FETI4IMatrix	matrix,
		FETI4IMPIInt	neighbours_size,
		FETI4IMPIInt*	neighbours,
		FETI4IInt		dirichlet_size,
		FETI4IInt*		dirichlet_indices,
		FETI4IReal*		dirichlet_values,
		FETI4IInt*		integer_options,
		FETI4IReal*		real_options
);


/// Solve an instance by the ESPRESO solver
/**
 * The method gets an instance and returns the solution to a pre-allocated array
 *
 * @param instance a data holder of an ESPRESO instance
 * @param rhs right-hand side
 * @param solution an array where is save the solution
 */
void FETI4ISolve(
		FETI4IInstance	instance,          // pointer to an instance
		FETI4IReal*		rhs,
		FETI4IReal*		solution           // solution
);


/*-----------------------------------------------------------------------------
 Functions for updating a created instance
------------------------------------------------------------------------------*/

void FETI4IUpdateStiffnessMatrix(
		FETI4IInstance	instance,
		FETI4IMatrix	stiffnessMatrix
);

void FETI4IUpdateRhs(
		FETI4IInstance	instance,
		FETI4IInt		size,
		FETI4IReal*		values
);

void FETI4IUpdateDirichlet(
		FETI4IInstance	instance,
		FETI4IInt		size,
		FETI4IInt*		indices,
		FETI4IReal*		values
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
		void*			ptr
);

#ifdef __cplusplus
}
#endif

/// List of allowed ESPRESO integer options.
/**
 * ESPRESO is controlled by passing parameters to FETI4ICreateInstance.
 * Integer parameters are in an array of size FETI4I_INTEGER_OPTIONS_SIZE.
 * See ESPRESO documentation for detailed description of each parameter.
 */
typedef enum {
	FETI4I_SUBDOMAINS,

	FETI4I_MAX_ITERATIONS,
	FETI4I_FETI_METHOD,
	FETI4I_PRECONDITIONER,
	FETI4I_CGSOLVER,

	FETI4I_INTEGER_OPTIONS_SIZE
} FETI4IIntegerOptions;

/// List of allowed ESPRESO real options.
/**
 * ESPRESO is controlled by passing parameters to FETI4ICreateInstance.
 * Real parameters are in an array of size FETI4I_REAL_OPTIONS_SIZE.
 * See ESPRESO documentation for detailed description of each parameter.
 */
typedef enum {
	FETI4I_PRECISION,

	FETI4I_REAL_OPTIONS_SIZE
} FETI4IRealOptions;

/// List of allowed matrix types
/**
 * ESPRESO matrix type allows some optimization of solver.
 */
typedef enum {
	FETI4I_REAL_SYMMETRIC_POSITIVE_DEFINITE,
	FETI4I_REAL_SYMMETRIC_INDEFINITE,
	FETI4I_REAL_UNSYMMETRIC,
} FETI4IMatrixType;

/// List of allowed element types
/**
 * ESPRESO needs to know a type of mesh elements.
 * The type corresponds with element dimension.
 */
typedef enum {
	FETI4I_ETYPE_POINT,
	FETI4I_ETYPE_LINE,
	FETI4I_ETYPE_PLANE,
	FETI4I_ETYPE_VOLUME,
} FETI4IElementType;

/// List of allowed DOFs orders
/**
 * The type corresponds with element dimension.
 */
typedef enum {
	FETI4I_ORDER_GROUPED_NODES,
	FETI4I_ORDER_GROUPED_DOFS
} FETI4IDOFsOrder;

#endif /* FETI4I_H_ */

