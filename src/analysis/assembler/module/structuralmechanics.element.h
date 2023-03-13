
#ifndef SRC_ANALYSIS_ASSEMBLER_MODULE_STRUCTURALMECHANICS_ELEMENT_H_
#define SRC_ANALYSIS_ASSEMBLER_MODULE_STRUCTURALMECHANICS_ELEMENT_H_

#include "math/simd/simd.h"

namespace espreso {

struct StructuralMechanicsGPC {
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

struct StructuralMechanicsElementType {
	enum: int {
		SYMMETRIC_PLANE              = 0,
		SYMMETRIC_PLANE_AXISYMMETRIC = 1,
		SYMMETRIC_VOLUME             = 2,
		FACE                         = 4,
		EDGE                         = 5,
		EDGE_AXISYMMETRIC            = 6,
		NODE                         = 7
	};
};

template <size_t nodes, size_t gps, size_t ndim, size_t edim, size_t etype> struct StructuralMechanicsDataDescriptor;

template <size_t nodes, size_t gps, size_t edim>
struct StructuralMechanicsDataDescriptor<nodes, gps, 2, edim, StructuralMechanicsElementType::SYMMETRIC_PLANE> {
	virtual ~StructuralMechanicsDataDescriptor() {}

	struct Element {
		struct {
			alignas(SIMD::size * sizeof(double)) SIMD thickness[gps];

			alignas(SIMD::size * sizeof(double)) SIMD youngModulus[gps];
			alignas(SIMD::size * sizeof(double)) SIMD poissonRatio[gps];
			alignas(SIMD::size * sizeof(double)) SIMD center      [gps][2]; // or rotation in the case of cartesion
			alignas(SIMD::size * sizeof(double)) SIMD elasticity  [gps][9];

			alignas(SIMD::size * sizeof(double)) SIMD density     [gps];
			alignas(SIMD::size * sizeof(double)) SIMD heatCapacity[gps];

			alignas(SIMD::size * sizeof(double)) SIMD acceleration   [gps][2];
			alignas(SIMD::size * sizeof(double)) SIMD angularVelocity[gps];
		} ecf;

		alignas(SIMD::size * sizeof(double)) SIMD coords[nodes][2];
		alignas(SIMD::size * sizeof(double)) SIMD gpcoords[gps][2];

		alignas(SIMD::size * sizeof(double)) double  w [gps];
		alignas(SIMD::size * sizeof(double)) double  N [gps][nodes];
		alignas(SIMD::size * sizeof(double)) double dN [gps][nodes][edim];

		alignas(SIMD::size * sizeof(double)) SIMD dND[gps][nodes][edim];
		alignas(SIMD::size * sizeof(double)) SIMD det[gps];

		alignas(SIMD::size * sizeof(double)) SIMD cossin    [gps][4];
		alignas(SIMD::size * sizeof(double)) SIMD elasticity[gps][9];

		Element()
		{
			for (size_t gp = 0; gp < gps; ++gp) {
				ecf.thickness[gp] = load1(1);
				ecf.density[gp] = load1(1);
				ecf.heatCapacity[gp] = load1(1);
			}
		}
	};

	virtual void simd(Element &element) =0;
	virtual void peel(Element &element, size_t size) { simd(element); }
};

template <size_t nodes, size_t gps, size_t edim>
struct StructuralMechanicsDataDescriptor<nodes, gps, 2, edim, StructuralMechanicsElementType::SYMMETRIC_PLANE_AXISYMMETRIC> {
	virtual ~StructuralMechanicsDataDescriptor() {}

	struct Element {
		struct {
			alignas(SIMD::size * sizeof(double)) SIMD thickness[gps];

			alignas(SIMD::size * sizeof(double)) SIMD youngModulus[gps];
			alignas(SIMD::size * sizeof(double)) SIMD poissonRatio[gps];
			alignas(SIMD::size * sizeof(double)) SIMD center      [gps][2]; // or rotation in the case of cartesion
			alignas(SIMD::size * sizeof(double)) SIMD elasticity  [gps][16];

			alignas(SIMD::size * sizeof(double)) SIMD density     [gps];
			alignas(SIMD::size * sizeof(double)) SIMD heatCapacity[gps];

			alignas(SIMD::size * sizeof(double)) SIMD acceleration   [gps][2];
			alignas(SIMD::size * sizeof(double)) SIMD angularVelocity[gps];
		} ecf;

		alignas(SIMD::size * sizeof(double)) SIMD coords[nodes][2];
		alignas(SIMD::size * sizeof(double)) SIMD gpcoords[gps][2];

		alignas(SIMD::size * sizeof(double)) double  w [gps];
		alignas(SIMD::size * sizeof(double)) double  N [gps][nodes];
		alignas(SIMD::size * sizeof(double)) double dN [gps][nodes][edim];

		alignas(SIMD::size * sizeof(double)) SIMD dND[gps][nodes][edim];
		alignas(SIMD::size * sizeof(double)) SIMD det[gps];

		alignas(SIMD::size * sizeof(double)) SIMD cossin    [gps][4];
		alignas(SIMD::size * sizeof(double)) SIMD elasticity[gps][16];

		Element()
		{
			for (size_t gp = 0; gp < gps; ++gp) {
				ecf.thickness[gp] = load1(1);
				ecf.density[gp] = load1(1);
				ecf.heatCapacity[gp] = load1(1);
			}
		}
	};

	virtual void simd(Element &element) =0;
	virtual void peel(Element &element, size_t size) { simd(element); }
};

template <size_t nodes, size_t gps, size_t edim>
struct StructuralMechanicsDataDescriptor<nodes, gps, 3, edim, StructuralMechanicsElementType::SYMMETRIC_VOLUME> {
	virtual ~StructuralMechanicsDataDescriptor() {}

	struct Element {
		struct {
			alignas(SIMD::size * sizeof(double)) SIMD youngModulus[gps][3];
			alignas(SIMD::size * sizeof(double)) SIMD poissonRatio[gps][3];
			alignas(SIMD::size * sizeof(double)) SIMD shearModulus[gps][3];
			alignas(SIMD::size * sizeof(double)) SIMD center      [gps][3]; // or rotation in the case of cartesion
			alignas(SIMD::size * sizeof(double)) SIMD elasticity  [gps][36];

			alignas(SIMD::size * sizeof(double)) SIMD density     [gps];
			alignas(SIMD::size * sizeof(double)) SIMD heatCapacity[gps];

			alignas(SIMD::size * sizeof(double)) SIMD acceleration   [gps][3];
			alignas(SIMD::size * sizeof(double)) SIMD angularVelocity[gps][3];
		} ecf;

		alignas(SIMD::size * sizeof(double)) SIMD coords[nodes][3];
		alignas(SIMD::size * sizeof(double)) SIMD gpcoords[gps][3];

		alignas(SIMD::size * sizeof(double)) double  w [gps];
		alignas(SIMD::size * sizeof(double)) double  N [gps][nodes];
		alignas(SIMD::size * sizeof(double)) double dN [gps][nodes][edim];

		alignas(SIMD::size * sizeof(double)) SIMD dND[gps][nodes][edim];
		alignas(SIMD::size * sizeof(double)) SIMD det[gps];

		alignas(SIMD::size * sizeof(double)) SIMD cossin    [gps][12];
		alignas(SIMD::size * sizeof(double)) SIMD elasticity[gps][36];

		Element()
		{
			for (size_t gp = 0; gp < gps; ++gp) {
				ecf.density[gp] = load1(1);
				ecf.heatCapacity[gp] = load1(1);
			}
		}
	};

	virtual void simd(Element &element) =0;
	virtual void peel(Element &element, size_t size) { simd(element); }
};

template <size_t nodes, size_t gps>
struct StructuralMechanicsDataDescriptor<nodes, gps, 3, 2, StructuralMechanicsElementType::FACE> {
	virtual ~StructuralMechanicsDataDescriptor() {}

	struct Element {
		struct {
			alignas(SIMD::size * sizeof(double)) SIMD normalPressure[gps];
		} ecf;

		alignas(SIMD::size * sizeof(double)) SIMD coords[nodes][3];
		alignas(SIMD::size * sizeof(double)) SIMD gpcoords[gps][3];

		alignas(SIMD::size * sizeof(double)) double  w [gps];
		alignas(SIMD::size * sizeof(double)) double  N [gps][nodes];
		alignas(SIMD::size * sizeof(double)) double dN [gps][nodes][2];

		alignas(SIMD::size * sizeof(double)) SIMD det[gps];
		alignas(SIMD::size * sizeof(double)) SIMD normal[gps][3];

		alignas(SIMD::size * sizeof(double)) SIMD displacement[nodes][3];

		Element()
		{

		}
	};

	virtual void simd(Element &element) =0;
	virtual void peel(Element &element, size_t size) { simd(element); }
};

template <size_t nodes, size_t ndim, size_t gps>
struct StructuralMechanicsDataDescriptor<nodes, gps, ndim, 1, StructuralMechanicsElementType::EDGE> {
	virtual ~StructuralMechanicsDataDescriptor() {}

	struct Element {
		struct {
			alignas(SIMD::size * sizeof(double)) SIMD normalPressure[gps];
		} ecf;

		alignas(SIMD::size * sizeof(double)) SIMD coords[nodes][ndim];
		alignas(SIMD::size * sizeof(double)) SIMD gpcoords[gps][ndim];

		alignas(SIMD::size * sizeof(double)) double  w [gps];
		alignas(SIMD::size * sizeof(double)) double  N [gps][nodes];
		alignas(SIMD::size * sizeof(double)) double dN [gps][nodes][1];

		alignas(SIMD::size * sizeof(double)) SIMD det[gps];
		alignas(SIMD::size * sizeof(double)) SIMD normal[gps][ndim];

		alignas(SIMD::size * sizeof(double)) SIMD displacement[nodes][ndim];

		Element()
		{

		}
	};

	virtual void simd(Element &element) =0;
	virtual void peel(Element &element, size_t size) { simd(element); }
};

template <size_t nodes, size_t ndim, size_t gps>
struct StructuralMechanicsDataDescriptor<nodes, gps, ndim, 1, StructuralMechanicsElementType::EDGE_AXISYMMETRIC> {
	virtual ~StructuralMechanicsDataDescriptor() {}

	struct Element {
		struct {
			alignas(SIMD::size * sizeof(double)) SIMD normalPressure[gps];
		} ecf;

		alignas(SIMD::size * sizeof(double)) SIMD coords[nodes][ndim];
		alignas(SIMD::size * sizeof(double)) SIMD gpcoords[gps][ndim];

		alignas(SIMD::size * sizeof(double)) double  w [gps];
		alignas(SIMD::size * sizeof(double)) double  N [gps][nodes];
		alignas(SIMD::size * sizeof(double)) double dN [gps][nodes][1];

		alignas(SIMD::size * sizeof(double)) SIMD det[gps];
		alignas(SIMD::size * sizeof(double)) SIMD normal[gps][ndim];

		alignas(SIMD::size * sizeof(double)) SIMD displacement[nodes][ndim];

		Element()
		{

		}
	};

	virtual void simd(Element &element) =0;
	virtual void peel(Element &element, size_t size) { simd(element); }
};

template <size_t ndim>
struct StructuralMechanicsDataDescriptor<1, 1, ndim, 0, StructuralMechanicsElementType::NODE> {
	virtual ~StructuralMechanicsDataDescriptor() {}

	struct Element {
		struct {
			alignas(SIMD::size * sizeof(double)) SIMD normalPressure[1];
		} ecf;

		alignas(SIMD::size * sizeof(double)) SIMD coords[1][ndim];
		alignas(SIMD::size * sizeof(double)) SIMD gpcoords[1][ndim];

		alignas(SIMD::size * sizeof(double)) SIMD displacement[1][ndim];

		Element()
		{

		}
	};

	virtual void simd(Element &element) =0;
	virtual void peel(Element &element, size_t size) { simd(element); }
};

}

#endif /* SRC_ANALYSIS_ASSEMBLER_MODULE_STRUCTURALMECHANICS_ELEMENT_H_ */
