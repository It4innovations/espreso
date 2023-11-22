
#ifndef SRC_ANALYSIS_ASSEMBLER_SUBKERNEL_SUBKERNEL_H_
#define SRC_ANALYSIS_ASSEMBLER_SUBKERNEL_SUBKERNEL_H_

#include "analysis/assembler/module/assembler.h"
#include "wrappers/simd/simd.h"

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

template <size_t size, size_t ndim> struct ElementCoordinatesSetter;
template <size_t size> struct ElementCoordinatesSetter<size, 2> {
	ElementCoordinatesSetter(SIMD data[][2], size_t i, size_t s, double &x, double &y, double &z)
	{
		x = data[i][0][s];
		y = data[i][1][s];
	}
};
template <size_t size> struct ElementCoordinatesSetter<size, 3> {
	ElementCoordinatesSetter(SIMD data[][3], size_t i, size_t s, double &x, double &y, double &z)
	{
		x = data[i][0][s];
		y = data[i][1][s];
		z = data[i][2][s];
	}
};

template <size_t nodes, size_t gps, size_t ndim> struct ElementCoordinates {
	alignas(SIMD::size * sizeof(double)) SIMD coords[nodes][ndim];
	alignas(SIMD::size * sizeof(double)) SIMD gpcoords[gps][ndim];

	void setCoordinatesNode(size_t n , size_t s, double &x, double &y, double &z) { ElementCoordinatesSetter<nodes, ndim>(coords  , n , s, x, y, z); }
	void setCoordinatesGP  (size_t gp, size_t s, double &x, double &y, double &z) { ElementCoordinatesSetter<gps  , ndim>(gpcoords, gp, s, x, y, z); }
};

template <size_t nodes, size_t gps> struct ElementTemperature {
	alignas(SIMD::size * sizeof(double)) SIMD temp[nodes];
	alignas(SIMD::size * sizeof(double)) SIMD gptemp[gps];

	void setTemperatureNode(size_t n , size_t s, double &temperature) { temperature = temp[ n][s]; }
	void setTemperatureGP  (size_t gp, size_t s, double &temperature) { temperature = temp[gp][s]; }
};

template <size_t nodes, size_t gps, size_t ndim> struct ElementDisplacement {
	alignas(SIMD::size * sizeof(double)) SIMD displacement[nodes][ndim];
	alignas(SIMD::size * sizeof(double)) SIMD smallStrainTensor[gps][ndim + (ndim * ndim - ndim) / 2];
};

template <size_t nodes, size_t gps, size_t edim> struct ElementIntegration {
	alignas(SIMD::size * sizeof(double)) double  w[gps];
	alignas(SIMD::size * sizeof(double)) double  N[gps][nodes];
	alignas(SIMD::size * sizeof(double)) double dN[gps][nodes][edim];

	alignas(SIMD::size * sizeof(double)) SIMD dND[gps][nodes][edim];
	alignas(SIMD::size * sizeof(double)) SIMD det[gps];
};

template <size_t nodes, size_t gps> struct ElementIntegration<nodes, gps, 0> {
	alignas(SIMD::size * sizeof(double)) double  w[gps];
	alignas(SIMD::size * sizeof(double)) double  N[gps][nodes];
};

template <size_t gps, size_t ndim> struct BondaryNormal {
	alignas(SIMD::size * sizeof(double)) SIMD normal[gps][ndim];
};

}



#endif /* SRC_ANALYSIS_ASSEMBLER_SUBKERNEL_SUBKERNEL_H_ */
