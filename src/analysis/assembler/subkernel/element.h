
#ifndef SRC_ANALYSIS_ASSEMBLER_SUBKERNEL_ELEMENT_H_
#define SRC_ANALYSIS_ASSEMBLER_SUBKERNEL_ELEMENT_H_

#include "math/simd/simd.h"

namespace espreso {

template <size_t gps, size_t nodes, size_t ndim> struct HeatTransferElement {
	alignas(SIMD::size * sizeof(double)) SIMD temp[nodes];
	alignas(SIMD::size * sizeof(double)) SIMD gptemp[gps];
	alignas(SIMD::size * sizeof(double)) SIMD coords[nodes][ndim];
	alignas(SIMD::size * sizeof(double)) SIMD gpcoords[gps][ndim];
};

template <size_t nodes, size_t gps, size_t edim> struct HeatTransferElementIntegration {
	alignas(SIMD::size * sizeof(double)) double  w[gps];
	alignas(SIMD::size * sizeof(double)) double  N[gps][nodes];
	alignas(SIMD::size * sizeof(double)) double dN[gps][nodes][edim];

	alignas(SIMD::size * sizeof(double)) SIMD dND[gps][nodes][edim];
	alignas(SIMD::size * sizeof(double)) SIMD det[gps];
};

}

#endif /* SRC_ANALYSIS_ASSEMBLER_SUBKERNEL_ELEMENT_H_ */
