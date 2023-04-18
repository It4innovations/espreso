
#ifndef SRC_ANALYSIS_ASSEMBLER_MODULE_HEATTRANSFER_ELEMENT_H_
#define SRC_ANALYSIS_ASSEMBLER_MODULE_HEATTRANSFER_ELEMENT_H_

#include "math/simd/simd.h"

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

struct HeatTransferElementType {
	enum: int {
		SYMMETRIC_ISOTROPIC  = 0,
		SYMMETRIC_GENERAL    = 1,
		ASYMMETRIC_ISOTROPIC = 2,
		ASYMMETRIC_GENERAL   = 3,
		FACE                 = 4,
		EDGE                 = 5,
		NODE                 = 6
	};
};

template <size_t gps, size_t ndim, enum ThermalConductivityConfiguration::MODEL> struct HeatTransferConductivity;

template <size_t gps, size_t ndim> struct HeatTransferConductivity<gps, ndim, ThermalConductivityConfiguration::MODEL::ISOTROPIC> {
	alignas(SIMD::size * sizeof(double)) SIMD conductivity[gps];
};

template <size_t gps, size_t ndim> struct HeatTransferConductivity<gps, ndim, ThermalConductivityConfiguration::MODEL::DIAGONAL> {
	alignas(SIMD::size * sizeof(double)) SIMD conductivity[gps][ndim];
};

template <size_t gps> struct HeatTransferConductivity<gps, 2, ThermalConductivityConfiguration::MODEL::SYMMETRIC> {
	alignas(SIMD::size * sizeof(double)) SIMD conductivity[gps][3];
};

template <size_t gps> struct HeatTransferConductivity<gps, 3, ThermalConductivityConfiguration::MODEL::SYMMETRIC> {
	alignas(SIMD::size * sizeof(double)) SIMD conductivity[gps][6];
};

template <size_t gps> struct HeatTransferConductivity<gps, 3, ThermalConductivityConfiguration::MODEL::ANISOTROPIC> {
	alignas(SIMD::size * sizeof(double)) SIMD conductivity[gps][ndim * ndim];
};

template <size_t gps, size_t ndim> struct HeatTransferElementParameters;

template <size_t gps> struct HeatTransferElementParameters<gps, 2> {
	alignas(SIMD::size * sizeof(double)) SIMD thickness   [gps];
	alignas(SIMD::size * sizeof(double)) SIMD density     [gps];
	alignas(SIMD::size * sizeof(double)) SIMD heatCapacity[gps];
	alignas(SIMD::size * sizeof(double)) SIMD heatSource  [gps];
	alignas(SIMD::size * sizeof(double)) SIMD advection   [gps][2];
};

template <size_t nodes, size_t gps, size_t ndim, size_t edim, size_t etype> struct HeatTransferDataDescriptor;

template <size_t nodes, size_t gps, size_t edim>
struct HeatTransferDataDescriptor<nodes, gps, 3, edim, HeatTransferElementType::SYMMETRIC_ISOTROPIC> {
	virtual ~HeatTransferDataDescriptor() {}

	struct Element {
		struct {
			alignas(SIMD::size * sizeof(double)) SIMD conductivity[gps];
			alignas(SIMD::size * sizeof(double)) SIMD density     [gps];
			alignas(SIMD::size * sizeof(double)) SIMD heatCapacity[gps];
			alignas(SIMD::size * sizeof(double)) SIMD heatSource  [gps];
		} ecf;

		alignas(SIMD::size * sizeof(double)) SIMD temp[nodes];
		alignas(SIMD::size * sizeof(double)) SIMD gptemp[gps];
		alignas(SIMD::size * sizeof(double)) SIMD coords[nodes][3];
		alignas(SIMD::size * sizeof(double)) SIMD gpcoords[gps][3];

		alignas(SIMD::size * sizeof(double)) double  w [gps];
		alignas(SIMD::size * sizeof(double)) double  N [gps][nodes];
		alignas(SIMD::size * sizeof(double)) double dN [gps][nodes][edim];

		alignas(SIMD::size * sizeof(double)) SIMD dND[gps][nodes][edim];
		alignas(SIMD::size * sizeof(double)) SIMD det[gps];

		alignas(SIMD::size * sizeof(double)) SIMD conductivity[gps];

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

template <size_t nodes, size_t gps, size_t edim>
struct HeatTransferDataDescriptor<nodes, gps, 2, edim, HeatTransferElementType::SYMMETRIC_ISOTROPIC> {
	virtual ~HeatTransferDataDescriptor() {}

	struct Element {
		struct {
			alignas(SIMD::size * sizeof(double)) SIMD thickness[gps];

			alignas(SIMD::size * sizeof(double)) SIMD conductivity[gps];
			alignas(SIMD::size * sizeof(double)) SIMD density     [gps];
			alignas(SIMD::size * sizeof(double)) SIMD heatCapacity[gps];
			alignas(SIMD::size * sizeof(double)) SIMD heatSource  [gps];
		} ecf;

		alignas(SIMD::size * sizeof(double)) SIMD temp[nodes];
		alignas(SIMD::size * sizeof(double)) SIMD gptemp[gps];
		alignas(SIMD::size * sizeof(double)) SIMD coords[nodes][2];
		alignas(SIMD::size * sizeof(double)) SIMD gpcoords[gps][2];

		alignas(SIMD::size * sizeof(double)) double  w [gps];
		alignas(SIMD::size * sizeof(double)) double  N [gps][nodes];
		alignas(SIMD::size * sizeof(double)) double dN [gps][nodes][edim];

		alignas(SIMD::size * sizeof(double)) SIMD dND[gps][nodes][edim];
		alignas(SIMD::size * sizeof(double)) SIMD det[gps];

		alignas(SIMD::size * sizeof(double)) SIMD conductivity[gps];

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
struct HeatTransferDataDescriptor<nodes, gps, 3, edim, HeatTransferElementType::SYMMETRIC_GENERAL> {
	virtual ~HeatTransferDataDescriptor() {}

	struct Element {
		struct {
			alignas(SIMD::size * sizeof(double)) SIMD conductivity[gps][6];
			alignas(SIMD::size * sizeof(double)) SIMD center      [gps][3]; // or rotation in the case of cartesion
			alignas(SIMD::size * sizeof(double)) SIMD density     [gps];
			alignas(SIMD::size * sizeof(double)) SIMD heatCapacity[gps];
			alignas(SIMD::size * sizeof(double)) SIMD heatSource  [gps];
		} ecf;

		alignas(SIMD::size * sizeof(double)) SIMD temp[nodes];
		alignas(SIMD::size * sizeof(double)) SIMD gptemp[gps];
		alignas(SIMD::size * sizeof(double)) SIMD coords[nodes][3];
		alignas(SIMD::size * sizeof(double)) SIMD gpcoords[gps][3];

		alignas(SIMD::size * sizeof(double)) double  w [gps];
		alignas(SIMD::size * sizeof(double)) double  N [gps][nodes];
		alignas(SIMD::size * sizeof(double)) double dN [gps][nodes][edim];

		alignas(SIMD::size * sizeof(double)) SIMD dND[gps][nodes][edim];
		alignas(SIMD::size * sizeof(double)) SIMD det[gps];

		alignas(SIMD::size * sizeof(double)) SIMD cossin      [gps][6];
		alignas(SIMD::size * sizeof(double)) SIMD conductivity[gps][6];

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

template <size_t nodes, size_t gps, size_t edim>
struct HeatTransferDataDescriptor<nodes, gps, 2, edim, HeatTransferElementType::SYMMETRIC_GENERAL> {
	virtual ~HeatTransferDataDescriptor() {}

	struct Element {
		struct {
			alignas(SIMD::size * sizeof(double)) SIMD thickness[gps];

			alignas(SIMD::size * sizeof(double)) SIMD conductivity[gps][3];
			alignas(SIMD::size * sizeof(double)) SIMD center      [gps][2]; // or rotation in the case of cartesion
			alignas(SIMD::size * sizeof(double)) SIMD density     [gps];
			alignas(SIMD::size * sizeof(double)) SIMD heatCapacity[gps];
			alignas(SIMD::size * sizeof(double)) SIMD heatSource  [gps];
		} ecf;

		alignas(SIMD::size * sizeof(double)) SIMD temp[nodes];
		alignas(SIMD::size * sizeof(double)) SIMD gptemp[gps];
		alignas(SIMD::size * sizeof(double)) SIMD coords[nodes][2];
		alignas(SIMD::size * sizeof(double)) SIMD gpcoords[gps][2];

		alignas(SIMD::size * sizeof(double)) double  w [gps];
		alignas(SIMD::size * sizeof(double)) double  N [gps][nodes];
		alignas(SIMD::size * sizeof(double)) double dN [gps][nodes][edim];

		alignas(SIMD::size * sizeof(double)) SIMD dND[gps][nodes][edim];
		alignas(SIMD::size * sizeof(double)) SIMD det[gps];

		alignas(SIMD::size * sizeof(double)) SIMD cossin      [gps][2];
		alignas(SIMD::size * sizeof(double)) SIMD conductivity[gps][3];

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
struct HeatTransferDataDescriptor<nodes, gps, 3, edim, HeatTransferElementType::ASYMMETRIC_ISOTROPIC> {
	virtual ~HeatTransferDataDescriptor() {}

	struct Element {
		struct {
			alignas(SIMD::size * sizeof(double)) SIMD conductivity[gps];
			alignas(SIMD::size * sizeof(double)) SIMD density     [gps];
			alignas(SIMD::size * sizeof(double)) SIMD heatCapacity[gps];
			alignas(SIMD::size * sizeof(double)) SIMD heatSource  [gps];
			alignas(SIMD::size * sizeof(double)) SIMD advection   [gps][3];
		} ecf;

		alignas(SIMD::size * sizeof(double)) SIMD temp[nodes];
		alignas(SIMD::size * sizeof(double)) SIMD gptemp[gps];
		alignas(SIMD::size * sizeof(double)) SIMD coords[nodes][3];
		alignas(SIMD::size * sizeof(double)) SIMD gpcoords[gps][3];

		alignas(SIMD::size * sizeof(double)) double  w [gps];
		alignas(SIMD::size * sizeof(double)) double  N [gps][nodes];
		alignas(SIMD::size * sizeof(double)) double dN [gps][nodes][edim];

		alignas(SIMD::size * sizeof(double)) SIMD dND[gps][nodes][edim];
		alignas(SIMD::size * sizeof(double)) SIMD det[gps];

		alignas(SIMD::size * sizeof(double)) SIMD conductivity[gps];

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

template <size_t nodes, size_t gps, size_t edim>
struct HeatTransferDataDescriptor<nodes, gps, 2, edim, HeatTransferElementType::ASYMMETRIC_ISOTROPIC> {
	virtual ~HeatTransferDataDescriptor() {}

	struct Element {
		struct {
			alignas(SIMD::size * sizeof(double)) SIMD thickness[gps];

			alignas(SIMD::size * sizeof(double)) SIMD conductivity[gps];
			alignas(SIMD::size * sizeof(double)) SIMD density     [gps];
			alignas(SIMD::size * sizeof(double)) SIMD heatCapacity[gps];
			alignas(SIMD::size * sizeof(double)) SIMD heatSource  [gps];
			alignas(SIMD::size * sizeof(double)) SIMD advection   [gps][2];
		} ecf;

		alignas(SIMD::size * sizeof(double)) SIMD temp[nodes];
		alignas(SIMD::size * sizeof(double)) SIMD gptemp[gps];
		alignas(SIMD::size * sizeof(double)) SIMD coords[nodes][2];
		alignas(SIMD::size * sizeof(double)) SIMD gpcoords[gps][2];

		alignas(SIMD::size * sizeof(double)) double  w [gps];
		alignas(SIMD::size * sizeof(double)) double  N [gps][nodes];
		alignas(SIMD::size * sizeof(double)) double dN [gps][nodes][edim];

		alignas(SIMD::size * sizeof(double)) SIMD dND[gps][nodes][edim];
		alignas(SIMD::size * sizeof(double)) SIMD det[gps];

		alignas(SIMD::size * sizeof(double)) SIMD conductivity[gps];

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
struct HeatTransferDataDescriptor<nodes, gps, 3, edim, HeatTransferElementType::ASYMMETRIC_GENERAL> {
	virtual ~HeatTransferDataDescriptor() {}

	struct Element {
		struct {
			alignas(SIMD::size * sizeof(double)) SIMD conductivity[gps][9];
			alignas(SIMD::size * sizeof(double)) SIMD center      [gps][3]; // or rotation in the case of cartesion
			alignas(SIMD::size * sizeof(double)) SIMD density     [gps];
			alignas(SIMD::size * sizeof(double)) SIMD heatCapacity[gps];
			alignas(SIMD::size * sizeof(double)) SIMD heatSource  [gps];
			alignas(SIMD::size * sizeof(double)) SIMD advection   [gps][3];
		} ecf;

		alignas(SIMD::size * sizeof(double)) SIMD temp[nodes];
		alignas(SIMD::size * sizeof(double)) SIMD gptemp[gps];
		alignas(SIMD::size * sizeof(double)) SIMD coords[nodes][3];
		alignas(SIMD::size * sizeof(double)) SIMD gpcoords[gps][3];

		alignas(SIMD::size * sizeof(double)) double  w [gps];
		alignas(SIMD::size * sizeof(double)) double  N [gps][nodes];
		alignas(SIMD::size * sizeof(double)) double dN [gps][nodes][edim];

		alignas(SIMD::size * sizeof(double)) SIMD dND[gps][nodes][edim];
		alignas(SIMD::size * sizeof(double)) SIMD det[gps];

		alignas(SIMD::size * sizeof(double)) SIMD cossin      [gps][6];
		alignas(SIMD::size * sizeof(double)) SIMD conductivity[gps][9];

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

template <size_t nodes, size_t gps, size_t edim>
struct HeatTransferDataDescriptor<nodes, gps, 2, edim, HeatTransferElementType::ASYMMETRIC_GENERAL> {
	virtual ~HeatTransferDataDescriptor() {}

	struct Element {
		struct {
			alignas(SIMD::size * sizeof(double)) SIMD thickness[gps];

			alignas(SIMD::size * sizeof(double)) SIMD conductivity[gps][4];
			alignas(SIMD::size * sizeof(double)) SIMD center      [gps][2]; // or rotation in the case of cartesion
			alignas(SIMD::size * sizeof(double)) SIMD density     [gps];
			alignas(SIMD::size * sizeof(double)) SIMD heatCapacity[gps];
			alignas(SIMD::size * sizeof(double)) SIMD heatSource  [gps];
			alignas(SIMD::size * sizeof(double)) SIMD advection   [gps][2];
		} ecf;

		alignas(SIMD::size * sizeof(double)) SIMD temp[nodes];
		alignas(SIMD::size * sizeof(double)) SIMD gptemp[gps];
		alignas(SIMD::size * sizeof(double)) SIMD coords[nodes][2];
		alignas(SIMD::size * sizeof(double)) SIMD gpcoords[gps][2];

		alignas(SIMD::size * sizeof(double)) double  w [gps];
		alignas(SIMD::size * sizeof(double)) double  N [gps][nodes];
		alignas(SIMD::size * sizeof(double)) double dN [gps][nodes][edim];

		alignas(SIMD::size * sizeof(double)) SIMD dND[gps][nodes][edim];
		alignas(SIMD::size * sizeof(double)) SIMD det[gps];

		alignas(SIMD::size * sizeof(double)) SIMD cossin      [gps][2];
		alignas(SIMD::size * sizeof(double)) SIMD conductivity[gps][4];

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

template <size_t nodes, size_t gps>
struct HeatTransferDataDescriptor<nodes, gps, 3, 2, HeatTransferElementType::FACE> {
	virtual ~HeatTransferDataDescriptor() {}

	struct Element {
		struct {
			alignas(SIMD::size * sizeof(double)) SIMD heatFlow[gps];
			alignas(SIMD::size * sizeof(double)) SIMD heatFlux[gps];
			alignas(SIMD::size * sizeof(double)) SIMD htc     [gps];
			alignas(SIMD::size * sizeof(double)) SIMD extTemp [gps];
		} ecf;

		alignas(SIMD::size * sizeof(double)) SIMD temp[nodes];
		alignas(SIMD::size * sizeof(double)) SIMD coords[nodes][3];
		alignas(SIMD::size * sizeof(double)) SIMD gpcoords[gps][3];

		alignas(SIMD::size * sizeof(double)) double  w [gps];
		alignas(SIMD::size * sizeof(double)) double  N [gps][nodes];
		alignas(SIMD::size * sizeof(double)) double dN [gps][nodes][2];

		alignas(SIMD::size * sizeof(double)) SIMD det[gps];

		Element()
		{

		}
	};

	virtual void simd(Element &element) =0;
	virtual void peel(Element &element, size_t size) { simd(element); }
};

template <size_t nodes, size_t gps>
struct HeatTransferDataDescriptor<nodes, gps, 3, 1, HeatTransferElementType::EDGE> {
	virtual ~HeatTransferDataDescriptor() {}

	struct Element {
		struct {
			alignas(SIMD::size * sizeof(double)) SIMD heatFlow[gps];
			alignas(SIMD::size * sizeof(double)) SIMD heatFlux[gps];
			alignas(SIMD::size * sizeof(double)) SIMD htc     [gps];
			alignas(SIMD::size * sizeof(double)) SIMD extTemp [gps];
		} ecf;

		alignas(SIMD::size * sizeof(double)) SIMD temp[nodes];
		alignas(SIMD::size * sizeof(double)) SIMD coords[nodes][3];
		alignas(SIMD::size * sizeof(double)) SIMD gpcoords[gps][3];

		alignas(SIMD::size * sizeof(double)) double  w [gps];
		alignas(SIMD::size * sizeof(double)) double  N [gps][nodes];
		alignas(SIMD::size * sizeof(double)) double dN [gps][nodes][1];

		alignas(SIMD::size * sizeof(double)) SIMD det[gps];

		Element()
		{

		}
	};

	virtual void simd(Element &element) =0;
	virtual void peel(Element &element, size_t size) { simd(element); }
};

template <size_t nodes, size_t gps>
struct HeatTransferDataDescriptor<nodes, gps, 2, 1, HeatTransferElementType::EDGE> {
	virtual ~HeatTransferDataDescriptor() {}

	struct Element {
		struct {
			alignas(SIMD::size * sizeof(double)) SIMD heatFlow[gps];
			alignas(SIMD::size * sizeof(double)) SIMD heatFlux[gps];
			alignas(SIMD::size * sizeof(double)) SIMD htc     [gps];
			alignas(SIMD::size * sizeof(double)) SIMD extTemp [gps];
		} ecf;

		alignas(SIMD::size * sizeof(double)) SIMD thickness[gps];
		alignas(SIMD::size * sizeof(double)) SIMD temp[nodes];
		alignas(SIMD::size * sizeof(double)) SIMD coords[nodes][2];
		alignas(SIMD::size * sizeof(double)) SIMD gpcoords[gps][2];

		alignas(SIMD::size * sizeof(double)) double  w [gps];
		alignas(SIMD::size * sizeof(double)) double  N [gps][nodes];
		alignas(SIMD::size * sizeof(double)) double dN [gps][nodes][1];

		alignas(SIMD::size * sizeof(double)) SIMD det[gps];

		Element()
		{

		}
	};

	virtual void simd(Element &element) =0;
	virtual void peel(Element &element, size_t size) { simd(element); }
};

template <>
struct HeatTransferDataDescriptor<1, 1, 3, 0, HeatTransferElementType::NODE> {
	virtual ~HeatTransferDataDescriptor() {}

	struct Element {
		struct {
			alignas(SIMD::size * sizeof(double)) SIMD heatFlow[1];
			alignas(SIMD::size * sizeof(double)) SIMD heatFlux[1];
			alignas(SIMD::size * sizeof(double)) SIMD htc     [1];
			alignas(SIMD::size * sizeof(double)) SIMD extTemp [1];
		} ecf;

		alignas(SIMD::size * sizeof(double)) SIMD temp[1];
		alignas(SIMD::size * sizeof(double)) SIMD coords[1][3];

		Element()
		{

		}
	};

	virtual void simd(Element &element) =0;
	virtual void peel(Element &element, size_t size) { simd(element); }
};

template <>
struct HeatTransferDataDescriptor<1, 1, 2, 0, HeatTransferElementType::NODE> {
	virtual ~HeatTransferDataDescriptor() {}

	struct Element {
		struct {
			alignas(SIMD::size * sizeof(double)) SIMD heatFlow[1];
			alignas(SIMD::size * sizeof(double)) SIMD heatFlux[1];
			alignas(SIMD::size * sizeof(double)) SIMD htc     [1];
			alignas(SIMD::size * sizeof(double)) SIMD extTemp [1];
		} ecf;

		alignas(SIMD::size * sizeof(double)) SIMD thickness[1];
		alignas(SIMD::size * sizeof(double)) SIMD temp[1];
		alignas(SIMD::size * sizeof(double)) SIMD coords[1][2];

		Element()
		{

		}
	};

	virtual void simd(Element &element) =0;
	virtual void peel(Element &element, size_t size) { simd(element); }
};

template <typename Element>
inline void setTemperature(Element &element, const size_t &n, const size_t &s, double &temperature)
{
	temperature = element.temp[n][s];
}

template <typename Element>
inline void setGPTemperature(Element &element, const size_t &gp, const size_t &s, double &temperature)
{
	temperature = element.gptemp[gp][s];
}

template <typename Element>
inline void setCoordinates(Element &element, const size_t &n, const size_t &s, double &coordinate_x, double &coordinate_y)
{
	coordinate_x = element.coords[n][0][s];
	coordinate_y = element.coords[n][1][s];
}

template <typename Element>
inline void setCoordinates(Element &element, const size_t &n, const size_t &s, double &coordinate_x, double &coordinate_y, double &coordinate_z)
{
	coordinate_x = element.coords[n][0][s];
	coordinate_y = element.coords[n][1][s];
	coordinate_z = element.coords[n][2][s];
}

template <typename Element>
inline void setGPCoordinates(Element &element, const size_t &gp, const size_t &s, double &coordinate_x, double &coordinate_y)
{
	coordinate_x = element.gpcoords[gp][0][s];
	coordinate_y = element.gpcoords[gp][1][s];
}

template <>
inline void setGPCoordinates<HeatTransferDataDescriptor<1, 1, 2, 0, HeatTransferElementType::NODE>::Element>(HeatTransferDataDescriptor<1, 1, 2, 0, HeatTransferElementType::NODE>::Element &element, const size_t &gp, const size_t &s, double &coordinate_x, double &coordinate_y)
{
	coordinate_x = element.coords[gp][0][s];
	coordinate_y = element.coords[gp][1][s];
}

template <typename Element>
inline void setGPCoordinates(Element &element, const size_t &gp, const size_t &s, double &coordinate_x, double &coordinate_y, double &coordinate_z)
{
	coordinate_x = element.gpcoords[gp][0][s];
	coordinate_y = element.gpcoords[gp][1][s];
	coordinate_z = element.gpcoords[gp][2][s];
}
template <>
inline void setGPCoordinates<HeatTransferDataDescriptor<1, 1, 3, 0, HeatTransferElementType::NODE>::Element>(HeatTransferDataDescriptor<1, 1, 3, 0, HeatTransferElementType::NODE>::Element &element, const size_t &gp, const size_t &s, double &coordinate_x, double &coordinate_y, double &coordinate_z)
{
	coordinate_x = element.coords[gp][0][s];
	coordinate_y = element.coords[gp][1][s];
	coordinate_z = element.coords[gp][2][s];
}

}



#endif /* SRC_ANALYSIS_ASSEMBLER_MODULE_HEATTRANSFER_ELEMENT_H_ */
