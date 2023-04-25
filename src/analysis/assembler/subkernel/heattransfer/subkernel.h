
#ifndef SRC_ANALYSIS_ASSEMBLER_SUBKERNEL_HEATTRANSFER_SUBKERNELS_H_
#define SRC_ANALYSIS_ASSEMBLER_SUBKERNEL_HEATTRANSFER_SUBKERNELS_H_

#include "analysis/assembler/subkernel/subkernel.h"
#include "config/ecf/material/thermalconductivity.h"

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

template <size_t gps, size_t ndim, enum ThermalConductivityConfiguration::MODEL> struct HeatTransferConductivity;

template <size_t gps, size_t ndim> struct HeatTransferConductivity<gps, ndim, ThermalConductivityConfiguration::MODEL::ISOTROPIC> {
	alignas(SIMD::size * sizeof(double)) SIMD conductivity[gps][1];
};

template <size_t gps, size_t ndim> struct HeatTransferConductivity<gps, ndim, ThermalConductivityConfiguration::MODEL::DIAGONAL> {
	alignas(SIMD::size * sizeof(double)) SIMD conductivity[gps][ndim];
};

template <size_t gps, size_t ndim> struct HeatTransferConductivity<gps, ndim, ThermalConductivityConfiguration::MODEL::SYMMETRIC> {
	alignas(SIMD::size * sizeof(double)) SIMD conductivity[gps][3 * (ndim - 1)];
};

template <size_t gps, size_t ndim> struct HeatTransferConductivity<gps, ndim, ThermalConductivityConfiguration::MODEL::ANISOTROPIC> {
	alignas(SIMD::size * sizeof(double)) SIMD conductivity[gps][ndim * ndim];
};

template <size_t gps, size_t ndim> struct HeatTransferRotation {
	alignas(SIMD::size * sizeof(double)) SIMD center[gps][ndim];
};

template <size_t gps, size_t ndim, enum ThermalConductivityConfiguration::MODEL ecfmodel> struct HeatTransferRotationMatrix {
	alignas(SIMD::size * sizeof(double)) SIMD cossin[gps][2 * ndim];
};

template <size_t gps, size_t ndim> struct HeatTransferRotationMatrix<gps, ndim, ThermalConductivityConfiguration::MODEL::ISOTROPIC> {

};

template <size_t gps, size_t ndim> struct HeatTransferElementParameters;

template <size_t gps> struct HeatTransferElementParameters<gps, 2> {
	alignas(SIMD::size * sizeof(double)) SIMD thickness   [gps];
	alignas(SIMD::size * sizeof(double)) SIMD density     [gps];
	alignas(SIMD::size * sizeof(double)) SIMD heatCapacity[gps];
	alignas(SIMD::size * sizeof(double)) SIMD heatSource  [gps];
	alignas(SIMD::size * sizeof(double)) SIMD advection   [gps][2];
};

template <size_t gps> struct HeatTransferElementParameters<gps, 3> {
	alignas(SIMD::size * sizeof(double)) SIMD density     [gps];
	alignas(SIMD::size * sizeof(double)) SIMD heatCapacity[gps];
	alignas(SIMD::size * sizeof(double)) SIMD heatSource  [gps];
	alignas(SIMD::size * sizeof(double)) SIMD advection   [gps][3];
};

template <size_t gps, size_t ndim> struct HeatTransferBoundaryParameters;

template <size_t gps> struct HeatTransferBoundaryParameters<gps, 2> {
	alignas(SIMD::size * sizeof(double)) SIMD thickness[gps];

	alignas(SIMD::size * sizeof(double)) SIMD heatFlow[gps];
	alignas(SIMD::size * sizeof(double)) SIMD heatFlux[gps];
	alignas(SIMD::size * sizeof(double)) SIMD htc     [gps];
	alignas(SIMD::size * sizeof(double)) SIMD extTemp [gps];
};

template <size_t gps> struct HeatTransferBoundaryParameters<gps, 3> {
	alignas(SIMD::size * sizeof(double)) SIMD heatFlow[gps];
	alignas(SIMD::size * sizeof(double)) SIMD heatFlux[gps];
	alignas(SIMD::size * sizeof(double)) SIMD htc     [gps];
	alignas(SIMD::size * sizeof(double)) SIMD extTemp [gps];
};

template <size_t nodes, size_t gps, size_t ndim, size_t edim, enum ThermalConductivityConfiguration::MODEL ecfmodel, enum ThermalConductivityConfiguration::MODEL model> struct HeatTransferElementDescriptor {
	virtual ~HeatTransferElementDescriptor() {}

	struct Element:
			ElementCoordinates<nodes, gps, ndim>,
			ElementTemperature<nodes, gps>,
			ElementIntegration<nodes, gps, edim>,
			HeatTransferConductivity<gps, ndim, model>,
			HeatTransferRotationMatrix<gps, ndim, ecfmodel>
	{
		struct ECF:
				HeatTransferConductivity<gps, ndim, ecfmodel>,
				HeatTransferRotation<gps, ndim>,
				HeatTransferElementParameters<gps, ndim>
		{

		} ecf;
	};

	virtual void simd(Element &element) =0;
};

template <size_t nodes, size_t gps, size_t ndim, size_t edim> struct HeatTransferBoundaryDescriptor {
	virtual ~HeatTransferBoundaryDescriptor() {}

	struct Element:
			ElementCoordinates<nodes, gps, ndim>,
			ElementTemperature<nodes, gps>,
			ElementIntegration<nodes, gps, edim>
	{
		struct ECF:
				HeatTransferBoundaryParameters<gps, ndim>
		{

		} ecf;
	};

	virtual void simd(Element &element) =0;
};

}

#endif /* SRC_ANALYSIS_ASSEMBLER_SUBKERNEL_HEATTRANSFER_SUBKERNELS_H_ */
