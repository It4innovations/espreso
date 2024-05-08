
#ifndef SRC_ANALYSIS_ASSEMBLER_HEATTRANSFER_ELEMENT_H_
#define SRC_ANALYSIS_ASSEMBLER_HEATTRANSFER_ELEMENT_H_

#include "analysis/assembler/general/element.h"

namespace espreso {

struct HeatTransferGPC {
    enum: int {
        POINT1    =  1,

        LINE2     =  2,
        LINE3     =  3,

        TRIANGLE3 =  6,
        SQUARE4   =  4,

        TRIANGLE6 =  6,
        SQUARE8   =  9,

        TETRA4    =  4,
        PYRAMID5  =  8,
        PRISMA6   =  9,
        HEXA8     =  8,

        TETRA10   = 15,
        PYRAMID13 = 14,
        PRISMA15  =  9,
        HEXA20    =  8,
    };
};

template <size_t nodes, size_t gps, size_t ndim, size_t edim> struct HeatTransferElement: public GeneralElement<nodes, gps, ndim, edim> {
    struct {
        alignas(SIMD::size * sizeof(double)) SIMD conductivity[ndim * ndim];
        alignas(SIMD::size * sizeof(double)) SIMD advection   [ndim];
        alignas(SIMD::size * sizeof(double)) SIMD heatSource;
        alignas(SIMD::size * sizeof(double)) SIMD density;
        alignas(SIMD::size * sizeof(double)) SIMD heatCapacity;
    } ecf;

    alignas(SIMD::size * sizeof(double)) SIMD conductivity[ndim * ndim];

    struct {
        alignas(SIMD::size * sizeof(double)) SIMD center[ndim];
        alignas(SIMD::size * sizeof(double)) SIMD cossin[ndim * 2];
    } rotation;

    alignas(SIMD::size * sizeof(double)) SIMD heat;
    alignas(SIMD::size * sizeof(double)) SIMD gradient[ndim];
    alignas(SIMD::size * sizeof(double)) SIMD flux    [ndim];

    alignas(SIMD::size * sizeof(double)) SIMD K[nodes * nodes];
    alignas(SIMD::size * sizeof(double)) SIMD M[nodes * nodes];
    alignas(SIMD::size * sizeof(double)) SIMD f[nodes];

    HeatTransferElement()
    {
        for (size_t i = 0; i < nodes * nodes; ++i) {
            K[i] = zeros();
            M[i] = zeros();
        }
    }
};

template <size_t nodes, size_t gps, size_t ndim, size_t edim> struct HeatTransferBoundary: public GeneralBoundary<nodes, gps, ndim, edim> {
    struct {
        alignas(SIMD::size * sizeof(double)) SIMD heatFlow;
        alignas(SIMD::size * sizeof(double)) SIMD heatFlux;
        alignas(SIMD::size * sizeof(double)) SIMD htc;
        alignas(SIMD::size * sizeof(double)) SIMD extTemp;
    } ecf;

    alignas(SIMD::size * sizeof(double)) SIMD f[nodes];
};

template <size_t ndim> struct HeatTransferNode: public GeneralNode<ndim> {

};

}

#endif /* SRC_ANALYSIS_ASSEMBLER_HEATTRANSFER_ELEMENT_H_ */
