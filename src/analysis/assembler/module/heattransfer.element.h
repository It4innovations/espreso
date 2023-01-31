
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
		SYMMETRIC_ISOTROPIC  = 0,
		SYMMETRIC_GENERAL    = 1,
		ASYMMETRIC_ISOTROPIC = 2,
		ASYMMETRIC_GENERAL   = 3,
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

template <size_t nodes, size_t gps, size_t ndim, size_t edim, size_t etype> struct HeatTransferDataDescriptor;

template <size_t nodes, size_t gps, size_t edim>
struct HeatTransferDataDescriptor<nodes, gps, 3, edim, HeatTransferElementType::SYMMETRIC_ISOTROPIC> {
	virtual ~HeatTransferDataDescriptor() {}

	struct Element {
		struct {
			alignas(SIMD::size * sizeof(double)) SIMD density     [gps];
			alignas(SIMD::size * sizeof(double)) SIMD heatCapacity[gps];
		} ecf;

		alignas(SIMD::size * sizeof(double)) SIMD coords[nodes][ndim];
		alignas(SIMD::size * sizeof(double)) SIMD gpcoords[gps][ndim];

		alignas(SIMD::size * sizeof(double)) SIMD  w [gps];
		alignas(SIMD::size * sizeof(double)) SIMD  N [gps][nodes];
		alignas(SIMD::size * sizeof(double)) SIMD dN [gps][nodes][edim];

		alignas(SIMD::size * sizeof(double)) SIMD dND[gps][nodes][edim];
		alignas(SIMD::size * sizeof(double)) SIMD det[gps];

		alignas(SIMD::size * sizeof(double)) SIMD conductivity[gps];

		Element()
		{
			for (size_t gp = 0; gp < gps; ++gp) {
				ecf.density[gp] = 1;
				ecf.heatCapacity[gp] = 1;
			}
		}
	};

	virtual void sisd(Element &element) =0;
	virtual void simd(Element &element) =0;
	virtual void peel(Element &element, size_t size) { simd(element); }
};

template <size_t nodes, size_t gps, size_t edim>
struct HeatTransferDataDescriptor<nodes, gps, 2, edim, HeatTransferElementType::SYMMETRIC_ISOTROPIC> {
	virtual ~HeatTransferDataDescriptor() {}

	struct Element {
		struct {
			alignas(SIMD::size * sizeof(double)) SIMD thickness[gps];

			alignas(SIMD::size * sizeof(double)) SIMD density     [gps];
			alignas(SIMD::size * sizeof(double)) SIMD heatCapacity[gps];
		} ecf;

		alignas(SIMD::size * sizeof(double)) SIMD coords[nodes][ndim];
		alignas(SIMD::size * sizeof(double)) SIMD gpcoords[gps][ndim];

		alignas(SIMD::size * sizeof(double)) SIMD  w [gps];
		alignas(SIMD::size * sizeof(double)) SIMD  N [gps][nodes];
		alignas(SIMD::size * sizeof(double)) SIMD dN [gps][nodes][edim];

		alignas(SIMD::size * sizeof(double)) SIMD dND[gps][nodes][edim];
		alignas(SIMD::size * sizeof(double)) SIMD det[gps];

		alignas(SIMD::size * sizeof(double)) SIMD conductivity[gps];

		Element()
		{
			for (size_t gp = 0; gp < gps; ++gp) {
				ecf.thickness[gp] = 1;
				ecf.density[gp] = 1;
				ecf.heatCapacity[gp] = 1;
			}
		}
	};

	virtual void sisd(Element &element) =0;
	virtual void simd(Element &element) =0;
	virtual void peel(Element &element, size_t size) { simd(element); }
};

template <size_t nodes, size_t gps, size_t edim>
struct HeatTransferDataDescriptor<nodes, gps, 3, edim, HeatTransferElementType::SYMMETRIC_GENERAL> {
	virtual ~HeatTransferDataDescriptor() {}

	struct Element {
		struct {
			alignas(SIMD::size * sizeof(double)) SIMD conductivity[gps][9];
			alignas(SIMD::size * sizeof(double)) SIMD angle       [gps][6];
			alignas(SIMD::size * sizeof(double)) SIMD density     [gps];
			alignas(SIMD::size * sizeof(double)) SIMD heatCapacity[gps];
		} ecf;

		alignas(SIMD::size * sizeof(double)) SIMD coords[nodes][ndim];
		alignas(SIMD::size * sizeof(double)) SIMD gpcoords[gps][ndim];

		alignas(SIMD::size * sizeof(double)) SIMD  w [gps];
		alignas(SIMD::size * sizeof(double)) SIMD  N [gps][nodes];
		alignas(SIMD::size * sizeof(double)) SIMD dN [gps][nodes][edim];

		alignas(SIMD::size * sizeof(double)) SIMD dND[gps][nodes][edim];
		alignas(SIMD::size * sizeof(double)) SIMD det[gps];

		alignas(SIMD::size * sizeof(double)) SIMD conductivity[gps][9];

		Element()
		{
			for (size_t gp = 0; gp < gps; ++gp) {
				std::fill(ecf.conductivity[gp], ecf.conductivity[gp] + 9, 1.);
				std::fill(ecf.angle[gp], ecf.angle[gp] + 6, 0.);
				ecf.density[gp] = 1;
				ecf.heatCapacity[gp] = 1;

				std::fill(conductivity[gp], conductivity[gp] + 9, 0.);
			}
		}
	};

	virtual void sisd(Element &element) =0;
	virtual void simd(Element &element) =0;
	virtual void peel(Element &element, size_t size) { simd(element); }
};

template <size_t nodes, size_t gps, size_t edim>
struct HeatTransferDataDescriptor<nodes, gps, 2, edim, HeatTransferElementType::SYMMETRIC_GENERAL> {
	virtual ~HeatTransferDataDescriptor() {}

	struct Element {
		struct {
			alignas(SIMD::size * sizeof(double)) SIMD thickness[gps];

			alignas(SIMD::size * sizeof(double)) SIMD conductivity[gps][4];
			alignas(SIMD::size * sizeof(double)) SIMD angle       [gps][2];
			alignas(SIMD::size * sizeof(double)) SIMD density     [gps];
			alignas(SIMD::size * sizeof(double)) SIMD heatCapacity[gps];
		} ecf;

		alignas(SIMD::size * sizeof(double)) SIMD coords[nodes][2];
		alignas(SIMD::size * sizeof(double)) SIMD gpcoords[gps][2];

		alignas(SIMD::size * sizeof(double)) SIMD  w [gps];
		alignas(SIMD::size * sizeof(double)) SIMD  N [gps][nodes];
		alignas(SIMD::size * sizeof(double)) SIMD dN [gps][nodes][edim];

		alignas(SIMD::size * sizeof(double)) SIMD dND[gps][nodes][edim];
		alignas(SIMD::size * sizeof(double)) SIMD det[gps];

		alignas(SIMD::size * sizeof(double)) SIMD conductivity[gps][4];

		Element()
		{
			for (size_t gp = 0; gp < gps; ++gp) {
				ecf.thickness[gp] = 1;

				std::fill(ecf.conductivity[gp], ecf.conductivity[gp] + 4, 0.);
				std::fill(ecf.angle[gp], ecf.angle[gp] + 2, 0.);
				ecf.density[gp] = 1;
				ecf.heatCapacity[gp] = 1;

				std::fill(conductivity[gp], conductivity[gp] + 4, 0.);
			}
		}
	};

	virtual void sisd(Element &element) =0;
	virtual void simd(Element &element) =0;
	virtual void peel(Element &element, size_t size) { simd(element); }
};

}



#endif /* SRC_ANALYSIS_ASSEMBLER_MODULE_HEATTRANSFER_ELEMENT_H_ */
