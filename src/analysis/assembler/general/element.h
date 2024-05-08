
#ifndef SRC_ANALYSIS_ASSEMBLER_GENERAL_ELEMENT_H_
#define SRC_ANALYSIS_ASSEMBLER_GENERAL_ELEMENT_H_

#include "wrappers/simd/simd.h"

namespace espreso {

template <size_t nodes, size_t gps, size_t ndim, size_t edim> struct GeneralElement {
    alignas(SIMD::size * sizeof(double)) double  w[gps];
    alignas(SIMD::size * sizeof(double)) double  N[gps][nodes];
    alignas(SIMD::size * sizeof(double)) double dN[gps][nodes][edim];

    alignas(SIMD::size * sizeof(double)) SIMD dND[nodes][edim];
    alignas(SIMD::size * sizeof(double)) SIMD det;
    alignas(SIMD::size * sizeof(double)) SIMD invJ[ndim * ndim];

    struct {
        alignas(SIMD::size * sizeof(double)) SIMD node[nodes];
        alignas(SIMD::size * sizeof(double)) SIMD gp;
    } thickness;

    struct {
        alignas(SIMD::size * sizeof(double)) SIMD node[nodes][ndim];
        alignas(SIMD::size * sizeof(double)) SIMD gp[ndim];
    } coords;

    struct {
        alignas(SIMD::size * sizeof(double)) SIMD initial[nodes];
        alignas(SIMD::size * sizeof(double)) SIMD node[nodes];
        alignas(SIMD::size * sizeof(double)) SIMD gp;
    } temperature;

    GeneralElement()
    {
        for (size_t n = 0; n < nodes; ++n) {
            thickness.node[n] = load1(1);
        }
        thickness.gp = load1(1);
    }
};

template <size_t nodes, size_t gps, size_t ndim, size_t edim> struct GeneralBoundary {
    alignas(SIMD::size * sizeof(double)) double  w[gps];
    alignas(SIMD::size * sizeof(double)) double  N[gps][nodes];
    alignas(SIMD::size * sizeof(double)) double dN[gps][nodes][edim];

    alignas(SIMD::size * sizeof(double)) SIMD dND[nodes][edim];
    alignas(SIMD::size * sizeof(double)) SIMD det;


    struct {
        alignas(SIMD::size * sizeof(double)) SIMD node[nodes];
        alignas(SIMD::size * sizeof(double)) SIMD gp;
    } thickness;

    struct {
        alignas(SIMD::size * sizeof(double)) SIMD node[nodes][ndim];
        alignas(SIMD::size * sizeof(double)) SIMD gp[ndim];
    } coords;

    struct {
        alignas(SIMD::size * sizeof(double)) SIMD node[nodes];
        alignas(SIMD::size * sizeof(double)) SIMD gp;
    } temperature;

    alignas(SIMD::size * sizeof(double)) SIMD normal[ndim];

    GeneralBoundary()
    {
        for (size_t n = 0; n < nodes; ++n) {
            thickness.node[n] = load1(1);
        }
        thickness.gp = load1(1);
    }
};

template <size_t ndim> struct GeneralNode {
    struct {
        alignas(SIMD::size * sizeof(double)) SIMD node[1][ndim];
    } coords;

    struct {
        alignas(SIMD::size * sizeof(double)) SIMD initial[1];
        alignas(SIMD::size * sizeof(double)) SIMD node[1];
    } temperature;
};

}


#endif /* SRC_ANALYSIS_ASSEMBLER_GENERAL_ELEMENT_H_ */
