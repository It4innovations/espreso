
#ifndef SRC_ANALYSIS_ASSEMBLER_SUBKERNEL_SUBKERNEL_H_
#define SRC_ANALYSIS_ASSEMBLER_SUBKERNEL_SUBKERNEL_H_

#include "analysis/assembler/module/assembler.h"
#include "math/simd/simd.h"

namespace espreso {

struct SubKernel {
	int isconst, isactive;
	Assembler::Action action;

	SubKernel(): isconst(1), isactive(0), action(Assembler::VOID) {}
	virtual ~SubKernel() {}

	void setActiveness(Assembler::Action action)
	{
		isactive = isactive && (this->action & action);
	}

	void setActiveness(Assembler::Action action, int guard)
	{
		isactive = guard && isactive && (this->action & action);
	}
};

template <size_t nodes, size_t gps, size_t ndim> struct ElementCoordinates {
	alignas(SIMD::size * sizeof(double)) SIMD coords[nodes][ndim];
	alignas(SIMD::size * sizeof(double)) SIMD gpcoords[gps][ndim];
};

template <size_t nodes, size_t gps> struct ElementTemperature {
	alignas(SIMD::size * sizeof(double)) SIMD temp[nodes];
	alignas(SIMD::size * sizeof(double)) SIMD gptemp[gps];
};

template <size_t nodes, size_t gps, size_t ndim> struct ElementDisplacement {
	alignas(SIMD::size * sizeof(double)) SIMD displacement[nodes][ndim];
	alignas(SIMD::size * sizeof(double)) SIMD dispTensor[gps][ndim * ndim / 2];
};

template <size_t nodes, size_t gps, size_t edim> struct ElementIntegration {
	alignas(SIMD::size * sizeof(double)) double  w[gps];
	alignas(SIMD::size * sizeof(double)) double  N[gps][nodes];
	alignas(SIMD::size * sizeof(double)) double dN[gps][nodes][edim];

	alignas(SIMD::size * sizeof(double)) SIMD dND[gps][nodes][edim];
	alignas(SIMD::size * sizeof(double)) SIMD det[gps];
};

template <size_t gps, size_t ndim> struct BondaryNormal {
	alignas(SIMD::size * sizeof(double)) SIMD normal[gps][ndim];
};

}



#endif /* SRC_ANALYSIS_ASSEMBLER_SUBKERNEL_SUBKERNEL_H_ */
