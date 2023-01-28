
#ifndef SRC_ANALYSIS_ASSEMBLER_MODULE_HEATTRANSFER_ELEMENT_H_
#define SRC_ANALYSIS_ASSEMBLER_MODULE_HEATTRANSFER_ELEMENT_H_

#include "math/simd/simd.h"

namespace espreso {

struct HeatTransferGPC {
	enum: int {
		POINT1    =  0,

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

struct HeatTransferElementType {
	enum: int {
		SYMMETRIC_ISOTROPIC = 0,
		SYMMETRIC_GENERAL   = 1,
	};
};

inline void printElementData(const char* name, const double *values, size_t size)
{
	printf("%s", name);
	for (size_t i = 0; i < size; ++i) {
		if (i % SIMD::size == 0) {
			printf("\n");
		}
		printf(" %+f", values[i]);
	}
	printf("\n");
}

template <int nodes, int gps, int ndim, int edim, int etype>
struct HeatTransferDataDescriptor {
	virtual ~HeatTransferDataDescriptor() {}

	struct Element {
		struct {
			alignas(SIMD::size * sizeof(double)) double conductivity[SIMD::size * gps * 9];
			alignas(SIMD::size * sizeof(double)) double angle       [SIMD::size * gps * 6];
			alignas(SIMD::size * sizeof(double)) double density     [SIMD::size * gps];
			alignas(SIMD::size * sizeof(double)) double heatCapacity[SIMD::size * gps];
		} ecf;

		alignas(SIMD::size * sizeof(double)) double coords[SIMD::size * nodes * ndim];
		alignas(SIMD::size * sizeof(double)) double gpcoords[SIMD::size * gps * ndim];

		alignas(SIMD::size * sizeof(double)) double  w [SIMD::size * gps];
		alignas(SIMD::size * sizeof(double)) double  N [SIMD::size * gps * nodes];
		alignas(SIMD::size * sizeof(double)) double dN [SIMD::size * gps * nodes * edim];

		alignas(SIMD::size * sizeof(double)) double dND[SIMD::size * gps * nodes * edim];
		alignas(SIMD::size * sizeof(double)) double det[SIMD::size * gps];

		alignas(SIMD::size * sizeof(double)) double conductivity[SIMD::size * gps * 9]; // 3 x 3
//
//		Element()
//		{
//			std::fill(ecf.conductivity, ecf.conductivity + SIMD::size * gps * 9, 1);
//			std::fill(ecf.angle, ecf.angle+ SIMD::size * gps * 6, 0);
//			std::fill(ecf.density, ecf.density+ SIMD::size * gps, 1);
//			std::fill(ecf.heatCapacity, ecf.heatCapacity+ SIMD::size * gps, 1);
//		}
	};

	virtual void sisd(Element &element) =0;
	virtual void simd(Element &element) =0;
	virtual void peel(Element &element, size_t size) { simd(element); }
//
//	static void print(const Element &element)
//	{
//		printf("HeatTransferDataDescriptor<%d, %d, %d, %d, %d>\n", nodes, gps, ndim, edim, etype);
//		printElementData("ecf.conductivity", element.ecf.conductivity, SIMD::size * gps * 4);
//		printElementData("ecf.angle", element.ecf.angle, SIMD::size * gps * 4);
//		printElementData("ecf.density", element.ecf.density, SIMD::size * gps);
//		printElementData("ecf.heatCapacity", element.ecf.heatCapacity, SIMD::size * gps);
//
//		printElementData("coords", element.coords, SIMD::size * nodes * ndim);
//
//		printElementData("w", element.w, SIMD::size * gps);
//		printElementData("N", element.N, SIMD::size * gps * nodes);
//		printElementData("dN", element.dN, SIMD::size * gps * nodes * edim);
//
//		printElementData("dND", element.dND, SIMD::size * gps * nodes * edim);
//		printElementData("det", element.det, SIMD::size * gps);
//
//		printElementData("conductivity", element.conductivity, SIMD::size * gps * 4);
//	}
};

template <int nodes, int gps, int edim, int etype>
struct HeatTransferDataDescriptor<nodes, gps, 2, edim, etype> {
	virtual ~HeatTransferDataDescriptor() {}

	struct Element {
		struct {
			alignas(SIMD::size * sizeof(double)) double thickness[SIMD::size * gps];

			alignas(SIMD::size * sizeof(double)) double conductivity[SIMD::size * gps * 4];
			alignas(SIMD::size * sizeof(double)) double angle       [SIMD::size * gps * 4];
			alignas(SIMD::size * sizeof(double)) double density     [SIMD::size * gps];
			alignas(SIMD::size * sizeof(double)) double heatCapacity[SIMD::size * gps];
		} ecf;

		alignas(SIMD::size * sizeof(double)) double coords[SIMD::size * nodes * 2];
		alignas(SIMD::size * sizeof(double)) double gpcoords[SIMD::size * gps * 2];

		alignas(SIMD::size * sizeof(double)) double  w [SIMD::size * gps];
		alignas(SIMD::size * sizeof(double)) double  N [SIMD::size * gps * nodes];
		alignas(SIMD::size * sizeof(double)) double dN [SIMD::size * gps * nodes * edim];

		alignas(SIMD::size * sizeof(double)) double dND[SIMD::size * gps * nodes * edim];
		alignas(SIMD::size * sizeof(double)) double det[SIMD::size * gps];

		alignas(SIMD::size * sizeof(double)) double conductivity[SIMD::size * gps * 4]; // 2 x 2
//
//		Element()
//		{
//			std::fill(ecf.thickness, ecf.thickness + SIMD::size * gps, 1);
//			std::fill(ecf.conductivity, ecf.conductivity + SIMD::size * gps * 4, 1);
//			std::fill(ecf.angle, ecf.angle+ SIMD::size * gps * 4, 0);
//			std::fill(ecf.density, ecf.density+ SIMD::size * gps, 1);
//			std::fill(ecf.heatCapacity, ecf.heatCapacity+ SIMD::size * gps, 1);
//		}
	};

	virtual void sisd(Element &element) =0;
	virtual void simd(Element &element) =0;
	virtual void peel(Element &element, size_t size) { simd(element); }
//
//	static void print(const Element &element)
//	{
//		printf("HeatTransferDataDescriptor<%d, %d, %d, %d, %d>\n", nodes, gps, ndim, edim, etype);
//		printElementData("ecf.thickness", element.ecf.thickness, SIMD::size * gps * 4);
//		printElementData("ecf.conductivity", element.ecf.conductivity, SIMD::size * gps * 4);
//		printElementData("ecf.angle", element.ecf.angle, SIMD::size * gps * 4);
//		printElementData("ecf.density", element.ecf.density, SIMD::size * gps);
//		printElementData("ecf.heatCapacity", element.ecf.heatCapacity, SIMD::size * gps);
//
//		printElementData("coords", element.coords, SIMD::size * nodes * 2);
//
//		printElementData("w", element.w, SIMD::size * gps);
//		printElementData("N", element.N, SIMD::size * gps * nodes);
//		printElementData("dN", element.dN, SIMD::size * gps * nodes * edim);
//
//		printElementData("dND", element.dND, SIMD::size * gps * nodes * edim);
//		printElementData("det", element.det, SIMD::size * gps);
//
//		printElementData("conductivity", element.conductivity, SIMD::size * gps * 4);
//	}
};

}



#endif /* SRC_ANALYSIS_ASSEMBLER_MODULE_HEATTRANSFER_ELEMENT_H_ */
