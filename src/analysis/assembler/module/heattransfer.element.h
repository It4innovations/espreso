
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

template <size_t nodes, size_t gps, size_t ndim, size_t edim, size_t etype>
struct HeatTransferDataDescriptor {
	virtual ~HeatTransferDataDescriptor() {}

	struct Element {
		struct {
			alignas(SIMD::size * sizeof(double)) double conductivity[SIMD::size * gps * 9];
			alignas(SIMD::size * sizeof(double)) double angle       [SIMD::size * gps * 6];
			alignas(SIMD::size * sizeof(double)) double density     [SIMD::size * gps];
			alignas(SIMD::size * sizeof(double)) double heatCapacity[SIMD::size * gps];
		} ecf;

		alignas(SIMD::size * sizeof(double)) SIMD coords[nodes][ndim];
		alignas(SIMD::size * sizeof(double)) double gpcoords[SIMD::size * gps * ndim];

		alignas(SIMD::size * sizeof(double)) SIMD  w [gps];
		alignas(SIMD::size * sizeof(double)) SIMD  N [gps][nodes];
		alignas(SIMD::size * sizeof(double)) SIMD dN [gps][nodes][edim];

		alignas(SIMD::size * sizeof(double)) SIMD dND[gps][nodes][edim];
		alignas(SIMD::size * sizeof(double)) SIMD det[gps];

		alignas(SIMD::size * sizeof(double)) double conductivity[SIMD::size * gps * 9]; // 3 x 3

		Element()
		{
			std::fill(ecf.conductivity, ecf.conductivity + SIMD::size * gps * 9, 1);
			std::fill(ecf.angle, ecf.angle+ SIMD::size * gps * 6, 0);
			std::fill(ecf.density, ecf.density+ SIMD::size * gps, 1);
			std::fill(ecf.heatCapacity, ecf.heatCapacity+ SIMD::size * gps, 1);
		}
	};

	virtual void sisd(Element &element) =0;
	virtual void simd(Element &element) =0;
	virtual void peel(Element &element, size_t size) { simd(element); }
};

template <size_t nodes, size_t gps, size_t edim, size_t etype>
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

		alignas(SIMD::size * sizeof(double)) SIMD coords[nodes][2];
		alignas(SIMD::size * sizeof(double)) double gpcoords[SIMD::size * gps * 2];

		alignas(SIMD::size * sizeof(double)) SIMD  w [gps];
		alignas(SIMD::size * sizeof(double)) SIMD  N [gps][nodes];
		alignas(SIMD::size * sizeof(double)) SIMD dN [gps][nodes][edim];

		alignas(SIMD::size * sizeof(double)) SIMD dND[gps][nodes][edim];
		alignas(SIMD::size * sizeof(double)) SIMD det[gps];

		alignas(SIMD::size * sizeof(double)) double conductivity[SIMD::size * gps * 4]; // 2 x 2

		Element()
		{
			std::fill(ecf.thickness, ecf.thickness + SIMD::size * gps, 1);
			std::fill(ecf.conductivity, ecf.conductivity + SIMD::size * gps * 4, 1);
			std::fill(ecf.angle, ecf.angle+ SIMD::size * gps * 4, 0);
			std::fill(ecf.density, ecf.density+ SIMD::size * gps, 1);
			std::fill(ecf.heatCapacity, ecf.heatCapacity+ SIMD::size * gps, 1);
		}
	};

	virtual void sisd(Element &element) =0;
	virtual void simd(Element &element) =0;
	virtual void peel(Element &element, size_t size) { simd(element); }
};

}



#endif /* SRC_ANALYSIS_ASSEMBLER_MODULE_HEATTRANSFER_ELEMENT_H_ */
