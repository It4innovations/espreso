
#ifndef SRC_ANALYSIS_ASSEMBLER_OPERATORS_FLUX_H_
#define SRC_ANALYSIS_ASSEMBLER_OPERATORS_FLUX_H_

#include <analysis/assembler/subkernel/operator.h>
#include "analysis/assembler/module/heattransfer.element.h"
#include "esinfo/meshinfo.h"
#include "mesh/store/elementstore.h"

namespace espreso {

struct TemperatureFluxBase: ActionOperator {
	const char* name() const { return "TemperatureFluxBase"; }

	TemperatureFluxBase(size_t interval, NamedData *flux)
	: flux(flux->data.data() + info::mesh->dimension * info::mesh->elements->eintervals[interval].begin)
	{
		isconst = false;
		action = Action::SOLUTION;
	}

	double* flux;

	template <size_t nodes, size_t gps, class Element>
	void isotropic2D(Element &element, size_t size)
	{
		double * __restrict__ out = flux;
		SIMD f0 = zeros(), f1 = zeros();
		for (size_t gp = 0; gp < gps; ++gp) {
			SIMD fgp0 = zeros(), fgp1 = zeros();
			for (size_t n = 0; n < nodes; ++n) {
				fgp0 = fgp0 + element.dND[gp][n][0] * element.temp[n];
				fgp1 = fgp1 + element.dND[gp][n][1] * element.temp[n];
			}
			f0 = f0 + fgp0 * element.conductivity[gp];
			f1 = f1 + fgp1 * element.conductivity[gp];
		}
		SIMD scale = load1(1. / gps);
		f0 = f0 * scale;
		f1 = f1 * scale;
		for (size_t s = 0; s <size; ++s) {
			out[2 * s + 0] = f0[s];
			out[2 * s + 1] = f1[s];
		}
	}

	template <size_t nodes, size_t gps, class Element>
	void isotropic3D(Element &element, size_t size)
	{
		double * __restrict__ out = flux;
		SIMD f0 = zeros(), f1 = zeros(), f2 = zeros();
		for (size_t gp = 0; gp < gps; ++gp) {
			SIMD fgp0 = zeros(), fgp1 = zeros(), fgp2 = zeros();
			for (size_t n = 0; n < nodes; ++n) {
				fgp0 = fgp0 + element.dND[gp][n][0] * element.temp[n];
				fgp1 = fgp1 + element.dND[gp][n][1] * element.temp[n];
				fgp2 = fgp2 + element.dND[gp][n][2] * element.temp[n];
			}
			f0 = f0 + fgp0 * element.conductivity[gp];
			f1 = f1 + fgp1 * element.conductivity[gp];
			f2 = f2 + fgp2 * element.conductivity[gp];
		}
		SIMD scale = load1(1. / gps);
		f0 = f0 * scale;
		f1 = f1 * scale;
		f2 = f2 * scale;
		for (size_t s = 0; s < size; ++s) {
			out[3 * s + 0] = f0[s];
			out[3 * s + 1] = f1[s];
			out[3 * s + 2] = f2[s];
		}
	}
};

template <size_t nodes, size_t gps, size_t ndim, size_t edim, size_t etype, class Physics> struct TemperatureFlux;

template <size_t nodes, size_t gps, size_t edim, class Physics>
struct TemperatureFlux<nodes, gps, 2, edim, HeatTransferElementType::SYMMETRIC_ISOTROPIC, Physics>: TemperatureFluxBase, Physics {
	using TemperatureFluxBase::TemperatureFluxBase;

	void move(int n)
	{
		flux += 2 * n;
	}

	void simd(typename Physics::Element &element)
	{
		peel(element, SIMD::size);
	}

	void peel(typename Physics::Element &element, size_t size)
	{
		isotropic2D<nodes, gps, typename Physics::Element>(element, size);
		move(size);
	}
};

template <size_t nodes, size_t gps, size_t edim, class Physics>
struct TemperatureFlux<nodes, gps, 3, edim, HeatTransferElementType::SYMMETRIC_ISOTROPIC, Physics>: TemperatureFluxBase, Physics {
	using TemperatureFluxBase::TemperatureFluxBase;

	void move(int n)
	{
		flux += 3 * n;
	}

	void simd(typename Physics::Element &element)
	{
		peel(element, SIMD::size);
	}

	void peel(typename Physics::Element &element, size_t size)
	{
		isotropic3D<nodes, gps, typename Physics::Element>(element, size);
		move(size);
	}
};

template <size_t nodes, size_t gps, size_t edim, class Physics>
struct TemperatureFlux<nodes, gps, 2, edim, HeatTransferElementType::ASYMMETRIC_ISOTROPIC, Physics>: TemperatureFluxBase, Physics {
	using TemperatureFluxBase::TemperatureFluxBase;

	void move(int n)
	{
		flux += 2 * n;
	}

	void simd(typename Physics::Element &element)
	{
		peel(element, SIMD::size);
	}

	void peel(typename Physics::Element &element, size_t size)
	{
		isotropic2D<nodes, gps, typename Physics::Element>(element, size);
		move(size);
	}
};

template <size_t nodes, size_t gps, size_t edim, class Physics>
struct TemperatureFlux<nodes, gps, 3, edim, HeatTransferElementType::ASYMMETRIC_ISOTROPIC, Physics>: TemperatureFluxBase, Physics {
	using TemperatureFluxBase::TemperatureFluxBase;

	void move(int n)
	{
		flux += 3 * n;
	}

	void simd(typename Physics::Element &element)
	{
		peel(element, SIMD::size);
	}

	void peel(typename Physics::Element &element, size_t size)
	{
		isotropic3D<nodes, gps, typename Physics::Element>(element, size);
		move(size);
	}
};

template <size_t nodes, size_t gps, size_t edim, class Physics>
struct TemperatureFlux<nodes, gps, 2, edim, HeatTransferElementType::SYMMETRIC_GENERAL, Physics>: TemperatureFluxBase, Physics {
	using TemperatureFluxBase::TemperatureFluxBase;

	void move(int n)
	{
		flux += 2 * n;
	}

	void simd(typename Physics::Element &element)
	{
		peel(element, SIMD::size);
	}

	void peel(typename Physics::Element &element, size_t size)
	{
		double * __restrict__ out = flux;
		SIMD f0 = zeros(), f1 = zeros();
		for (size_t gp = 0; gp < gps; ++gp) {
			SIMD fgp0 = zeros(), fgp1 = zeros();
			for (size_t n = 0; n < nodes; ++n) {
				fgp0 = fgp0 + element.dND[gp][n][0] * element.temp[n];
				fgp1 = fgp1 + element.dND[gp][n][1] * element.temp[n];
			}
			f0 = f0 + fgp0 * element.conductivity[gp][0] + fgp1 * element.conductivity[gp][1];
			f1 = f1 + fgp0 * element.conductivity[gp][1] + fgp1 * element.conductivity[gp][2];
		}
		SIMD scale = load1(1. / gps);
		f0 = f0 * scale;
		f1 = f1 * scale;
		for (size_t s = 0; s <size; ++s) {
			out[2 * s + 0] = f0[s];
			out[2 * s + 1] = f1[s];
		}
		move(size);
	}
};

template <size_t nodes, size_t gps, size_t edim, class Physics>
struct TemperatureFlux<nodes, gps, 2, edim, HeatTransferElementType::ASYMMETRIC_GENERAL, Physics>: TemperatureFluxBase, Physics {
	using TemperatureFluxBase::TemperatureFluxBase;

	void move(int n)
	{
		flux += 2 * n;
	}

	void simd(typename Physics::Element &element)
	{
		peel(element, SIMD::size);
	}

	void peel(typename Physics::Element &element, size_t size)
	{
		double * __restrict__ out = flux;
		SIMD f0 = zeros(), f1 = zeros();
		for (size_t gp = 0; gp < gps; ++gp) {
			SIMD fgp0 = zeros(), fgp1 = zeros();
			for (size_t n = 0; n < nodes; ++n) {
				fgp0 = fgp0 + element.dND[gp][n][0] * element.temp[n];
				fgp1 = fgp1 + element.dND[gp][n][1] * element.temp[n];
			}
			f0 = f0 + fgp0 * element.conductivity[gp][0] + fgp1 * element.conductivity[gp][1];
			f1 = f1 + fgp0 * element.conductivity[gp][2] + fgp1 * element.conductivity[gp][3];
		}
		SIMD scale = load1(1. / gps);
		f0 = f0 * scale;
		f1 = f1 * scale;
		for (size_t s = 0; s <size; ++s) {
			out[2 * s + 0] = f0[s];
			out[2 * s + 1] = f1[s];
		}
		move(size);
	}
};

template <size_t nodes, size_t gps, size_t edim, class Physics>
struct TemperatureFlux<nodes, gps, 3, edim, HeatTransferElementType::SYMMETRIC_GENERAL, Physics>: TemperatureFluxBase, Physics {
	using TemperatureFluxBase::TemperatureFluxBase;

	void move(int n)
	{
		flux += 3 * n;
	}

	void simd(typename Physics::Element &element)
	{
		peel(element, SIMD::size);
	}

	void peel(typename Physics::Element &element, size_t size)
	{
		double * __restrict__ out = flux;
		SIMD f0 = zeros(), f1 = zeros(), f2 = zeros();
		for (size_t gp = 0; gp < gps; ++gp) {
			SIMD fgp0 = zeros(), fgp1 = zeros(), fgp2 = zeros();
			for (size_t n = 0; n < nodes; ++n) {
				fgp0 = fgp0 + element.dND[gp][n][0] * element.temp[n];
				fgp1 = fgp1 + element.dND[gp][n][1] * element.temp[n];
				fgp2 = fgp2 + element.dND[gp][n][2] * element.temp[n];
			}
			f0 = f0 + fgp0 * element.conductivity[gp][0] + fgp1 * element.conductivity[gp][1] + fgp2 * element.conductivity[gp][2];
			f1 = f1 + fgp0 * element.conductivity[gp][1] + fgp1 * element.conductivity[gp][3] + fgp2 * element.conductivity[gp][4];
			f2 = f2 + fgp0 * element.conductivity[gp][2] + fgp1 * element.conductivity[gp][4] + fgp2 * element.conductivity[gp][5];
		}
		SIMD scale = load1(1. / gps);
		f0 = f0 * scale;
		f1 = f1 * scale;
		f2 = f2 * scale;
		for (size_t s = 0; s < size; ++s) {
			out[3 * s + 0] = f0[s];
			out[3 * s + 1] = f1[s];
			out[3 * s + 2] = f2[s];
		}
		move(size);
	}
};

template <size_t nodes, size_t gps, size_t edim, class Physics>
struct TemperatureFlux<nodes, gps, 3, edim, HeatTransferElementType::ASYMMETRIC_GENERAL, Physics>: TemperatureFluxBase, Physics {
	using TemperatureFluxBase::TemperatureFluxBase;

	void move(int n)
	{
		flux += 3 * n;
	}

	void simd(typename Physics::Element &element)
	{
		peel(element, SIMD::size);
	}

	void peel(typename Physics::Element &element, size_t size)
	{
		double * __restrict__ out = flux;
		SIMD f0 = zeros(), f1 = zeros(), f2 = zeros();
		for (size_t gp = 0; gp < gps; ++gp) {
			SIMD fgp0 = zeros(), fgp1 = zeros(), fgp2 = zeros();
			for (size_t n = 0; n < nodes; ++n) {
				fgp0 = fgp0 + element.dND[gp][n][0] * element.temp[n];
				fgp1 = fgp1 + element.dND[gp][n][1] * element.temp[n];
				fgp2 = fgp2 + element.dND[gp][n][2] * element.temp[n];
			}
			f0 = f0 + fgp0 * element.conductivity[gp][0] + fgp1 * element.conductivity[gp][1] + fgp2 * element.conductivity[gp][2];
			f1 = f1 + fgp0 * element.conductivity[gp][3] + fgp1 * element.conductivity[gp][4] + fgp2 * element.conductivity[gp][5];
			f2 = f2 + fgp0 * element.conductivity[gp][6] + fgp1 * element.conductivity[gp][7] + fgp2 * element.conductivity[gp][8];
		}
		SIMD scale = load1(1. / gps);
		f0 = f0 * scale;
		f1 = f1 * scale;
		f2 = f2 * scale;
		for (size_t s = 0; s < size; ++s) {
			out[3 * s + 0] = f0[s];
			out[3 * s + 1] = f1[s];
			out[3 * s + 2] = f2[s];
		}
		move(size);
	}
};

}

#endif /* SRC_ANALYSIS_ASSEMBLER_OPERATORS_FLUX_H_ */
