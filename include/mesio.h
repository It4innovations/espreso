
#ifndef _MESIO_H_
#define _MESIO_H_

#include "mpi.h"

/*-----------------------------------------------------------------------------
 Set data-types used in MESIO

 Possible values for MESIO_INT_WIDTH:
    32: use 32 bit signed integer
    64: use 64 bit signed integer

 Possible values for MESIO_REAL_WIDTH:
    64: ESPRESO supports only 64 bit real values
------------------------------------------------------------------------------*/
#ifndef MESIO_INT_WIDTH
#define MESIO_INT_WIDTH 32
#endif

#ifndef MESIO_REAL_WIDTH
#define MESIO_REAL_WIDTH 64
#endif

#if MESIO_INT_WIDTH == 32
    typedef int MESIOInt;
#elif MESIO_INT_WIDTH == 64
    typedef long MESIOInt;
#else
#error "Incorrect user-supplied value of MESIO_INT_WIDTH"
#endif

#if MESIO_REAL_WIDTH == 64
    typedef double MESIOReal;
#else
#error "Incorrect user-supplied value of MESIO_REAL_WIDTH"
#endif


typedef enum {
    MESIO_ANSYS,
    MESIO_ENSIGHT,
    MESIO_VTK_LEGACY,
    MESIO_XDMF,
} MESIOFormat;

typedef enum {
    MESIO_NONE,
    MESIO_METIS,
    MESIO_PARMETIS,
    MESIO_PTSCOTCH,
    MESIO_HILBERT_CURVE
} MESIODecomposer;

typedef enum {
    POINT1, // 0

    // without mid-points
    LINE2, // 1
    TRIANGLE3, // 2
    SQUARE4, // 3
    TETRA4, // 4
    PYRAMID5, // 5
    PRISMA6, // 6
    HEXA8, // 7

    // with mid-points
    LINE3, // 8
    TRIANGLE6, // 9
    SQUARE8, // 10
    TETRA10, // 11
    PYRAMID13, // 12
    PRISMA15, // 13
    HEXA20, // 14

    // element with unknown number of nodes
    NOT_SUPPORTED,

    SIZE
} MESIOElementType;

#ifdef __cplusplus
extern "C" {
#endif

/*-----------------------------------------------------------------------------
 Definitions of internal structures used in MESIO
------------------------------------------------------------------------------*/
typedef struct MESIOData* MESIO;

/*-----------------------------------------------------------------------------
 Functions for manipulating with MESIO internal structures
------------------------------------------------------------------------------*/

/// Initialize the MESIO library
/**
 * ESPRESO needs to be initialized before calling any other function.
 * At the end of the program MESIOFinalize should be called.
 *
 * @param comm an MPI communication
 * @param verbosity a verbosity level
 */
void MESIOInit(
    MPI_Comm        comm,
    int             verbosity
);

/// Finalize the MESIO library
/**
 * This function destroy all internal parameters.
 */
void MESIOFinalize();

void MESIOLoad(
    MESIO*          mesio,
    MESIOFormat     format,
    const char*     path,
    MESIODecomposer decomposer,
    int             domains
);

void MESIONodes(
    MESIO           mesio,
    MESIOInt*       nhalo,
    MESIOInt*       offset,
    MESIOInt*       size,
    MESIOInt*       totalSize,
    MESIOInt**      ids,
    MESIOInt**      position,
    MESIOReal**     coordinates
);

void MESIONodesRanks(
    MESIO           mesio,
    MESIOInt**      rankDistribution,
    int**           rankData
);

void MESIONodesDomains(
    MESIO           mesio,
    MESIOInt**      domainDistribution,
    MESIOInt**      domainData
);

void MESIOElements(
    MESIO           mesio,
    MESIOInt*       offset,
    MESIOInt*       size,
    MESIOInt*       totalSize,
    MESIOInt**      type,
    MESIOInt**      enodesDistribution,
    MESIOInt**      enodesData
);

void MESIOElementsDomains(
    MESIO           mesio,
	MESIOInt**      domains
);

void MESIOElementsMaterials(
    MESIO           mesio,
    int**           material
);

void MESIOElementsBodies(
    MESIO           mesio,
    MESIOInt*       bodies,
    int**           body
);

void MESIOElementsNeighbors(
    MESIO           mesio,
    MESIOInt**      neighborDistribution,
    MESIOInt**      neighborData
);

void MESIOElementsCounters(
    MESIO           mesio,
    MESIOInt        etype,
    MESIOInt*       offset,
    MESIOInt*       totalSize
);

int MESIOElementsRegions(
    MESIO           mesio
);

void MESIOElementsRegion(
    MESIO           mesio,
    MESIOInt        region,
    const char**    name,
    MESIOInt*       size,
    MESIOInt**      elements
);

void MESIOElementsRegionNodes(
    MESIO           mesio,
    MESIOInt        region,
    MESIOInt*       nhalo,
    MESIOInt*       offset,
    MESIOInt*       size,
    MESIOInt*       totalSize,
    MESIOInt**      nodes,
	MESIOInt**      position
);

void MESIOElementsRegionCounters(
    MESIO           mesio,
    MESIOInt        region,
    MESIOInt        etype,
    MESIOInt*       offset,
    MESIOInt*       totalSize
);

int MESIOBoundaryRegions(
    MESIO           mesio
);

void MESIOBoundaryRegion(
    MESIO           mesio,
    MESIOInt        region,
    const char**    name,
    MESIOInt*       dimension,
    MESIOInt*       size,
    MESIOInt**      type,
    MESIOInt**      parent,
    MESIOInt**      elementDistribution,
    MESIOInt**      elementData
);

void MESIOBoundaryRegionNodes(
    MESIO           mesio,
    MESIOInt        region,
    MESIOInt*       nhalo,
    MESIOInt*       offset,
    MESIOInt*       size,
    MESIOInt*       totalSize,
    MESIOInt**      nodes,
    MESIOInt**      position
);

void MESIOBoundaryRegionCounters(
    MESIO           mesio,
    MESIOInt        region,
    MESIOInt        etype,
    MESIOInt*       offset,
    MESIOInt*       totalSize
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
void MESIODestroy(
    void*        ptr
);

#ifdef __cplusplus
}
#endif

#endif /* _MESIO_H_ */
